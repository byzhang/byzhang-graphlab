/**
 * Copyright (c) 2009 Carnegie Mellon University.
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */


/**
 * \file
 *
 * This file contains the template for the redis based graphlab graph
 * data-structure.
 *
 */

#ifndef GRAPHLAB_REDIS_GRAPH_HPP
#define GRAPHLAB_REDIS_GRAPH_HPP

#include <cmath>

#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <algorithm>
#include <functional>
#include <fstream>

#include <boost/bind.hpp>
#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#include <boost/coerce.hpp>
#include <boost/unordered_set.hpp>

#include <hiredispp.h>

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/serialization/serialization-includes.hpp>
#include <graphlab/util/random.hpp>

#include <graphlab/macros_def.hpp>

using namespace boost;
using namespace hiredispp;

namespace graphlab {
  // CLASS GRAPH ==============================================================>
  /**
     \brief The GraphLab primary Graph container templatized over the
     vertex and edge types.

     Every vertex and edge in the graph is assigned a unique integer
     ID.  The type of the vertex id is
     <code>graphlab::vertex_id_type</code> and the type of the edge id
     is <code>graphlab::edge_id_type</code>. Both
     <code>vertex_id_type</code> and <code>edge_id_type</code> are
     currently defined as <code>uint32_t</code>.  While this limits
     the graphs to 4 billion vertices it also helps reduce the storage
     overhead. We encourage users to use the
     <code>vertex_id_type</code> and <code>edge_id_type</code> types
     as they may change in larger distributed systems.

     <h2> Graph Creation </h2>

     Vertices and edges are added using the graph::add_vertex()
     and graph::add_edge() member functions:


     \code
     vertex_id_type graph::add_vertex(const VertexData& vdata = VertexData())
     edge_id_type graph::add_edge(vertex_id_type source, vertex_id_type target,
                                  const EdgeData& edata = EdgeData())
     \endcode

     The functions return the id's of the added vertex and edge
     respectively.  An edge can only be added if both the source and
     target vertex id's are already in the graph. Duplicate edges are not
     supported and may result in undefined behavior.
  */
  template<typename VertexData, typename EdgeData>
  class redis_graph {
  public:

    /// The type of a vertex is a simple size_t
    typedef uint32_t vertex_id_type;

    /// The type of an edge id
    typedef uint32_t edge_id_type;

    /// Type for vertex colors
    typedef vertex_id_type vertex_color_type;

    /** The type of the edge list */
    typedef vector<edge_id_type> edge_list_type;

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData edge_data_type;

  public:
    struct redis_server {
      const char *host;
      int port;
      int db;
      redis_server(const char *host_, int port_=6379, int db_=0) {
        host = host_;
        port = port_;
        db = db_;
      }
    }

    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    redis_graph(const redis_server& vertex_server,
                const redis_server& edge_server,
                const redis_server& in_edge_server,
                const redis_server& out_edge_server,
                const redis_server& color_server)
      : vertices(new Redis(vertex_server.host, vertex_server.port))
      , vcolors(new Redis(color_server.host, color_server.port))
      , edges(new Redis(edge_server.host, edge_server.port))
      , in_edges(new Redis(in_edge_server.host, in_edge_server.port))
      , out_edges(new Redis(out_edge_server.host, out_edge_server.port))
    {
      // TODO: check/handle the same server address.
      vertices->select(vertex_server.db);
      vcolors->select(color_server.db);
      edges->select(edge_server.db);
      edges->select(edge_server.db);
      in_edges->select(in_edge_server.db);
      out_edges->select(out_edge_server.db);
    }

    virtual ~redis_graph() {
      delete vertices;
      delete vcolors;
      delete edges;
      delete in_edges;
      delete out_edges;
    }

    // METHODS =================================================================>

    /**
     * \brief Resets the graph state.
     */
    void clear() {
      vertices->flush();
      vcolors->flush();
      edges->flush();
      in_edges->flush();
      out_edges->flush();
      ++changeid;
    }

    /** \brief Get the number of vertices */
    size_t num_vertices() const {
      return vertices->size();
    } // end of num vertices

    /** \brief Get the number of vertices local to this machine */
    size_t local_vertices() const {
      return vertices->size();
    } // end of num vertices

    /** \brief Get the number of edges */
    size_t num_edges() const {
      return edges->size();
    } // end of num edges


    /** \brief Get the number of in edges of a particular vertex */
    size_t num_in_neighbors(vertex_id_type v) const {
      return in_edges->hlen(v);
    } // end of num vertices

    /** \brief Get the number of out edges of a particular vertex */
    size_t num_out_neighbors(vertex_id_type v) const  {
      return out_edges->hlen(v);
    } // end of num vertices

    /** \brief Finds an edge.
        The value of the first element of the pair will be true if an
        edge from src to target is found and false otherwise. If the
        edge is found, the edge ID is returned in the second element of the pair. */
    std::pair<bool, edge_id_type>
    find(vertex_id_type source, vertex_id_type target) const {
      auto result = out_edges->hget(source, target);
      if (result == RedisConst<char>::Nil) {
        return make_pair(false, -1);
      }
      return make_pair(true, HexString2Int<edge_id_type>(result));
    } // end of find


    /** \brief A less safe version of find.
        Returns the edge_id of an edge from src to target exists.
        Assertion failure otherwise. */
    edge_id_type edge_id(vertex_id_type source, vertex_id_type target) const {
      auto res = find(source, target);
      // The edge must exist
      ASSERT_TRUE(res.first);
      DCHECK_LT(res.second, edges->size());
      return res.second;
    } // end of edge_id


    /** \brief Returns the edge ID of the edge going in the opposite direction.
        Assertion failure if such an edge is not found.  */
    edge_id_type rev_edge_id(edge_id_type eid) const {
      DCHECK_LT(eid, edges->size());
      Redis::Reply reply = edges->hmget(eid, 'S', 'T');
      string source = reply[0];
      string target = Reply[1];
      return edge_id(HexString2Int<vertext_id_type>(target), HexString2Int(source));
    } // end of rev_edge_id

    /**
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     */
    vertex_id_type add_vertex(const VertexData& vdata = VertexData() ) {
      vertex_id_type id = 0;
      while (true) {
        id = vertices->size();
        bool ok = vertices->setnx(Int2HexString(id), serialize_to_string(vdata));
        if (ok) {
          break;
        }
      }
      ASSERT_EQ(1, vcolors->setnx(Int2HexString(id), 0));
      return id;
    } // End of add vertex;

    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    edge_id_type add_edge(vertex_id_type source, vertex_id_type target,
                          const EdgeData& edata = EdgeData()) {
      auto num_vertices = vertices->size();
      if (source >= num_vertices || target >= num_vertices) {
        logstream(LOG_FATAL)
          << "Attempting add_edge (" << source
          << " -> " << target
          << ") when there are only " << vertices->size()
          << " vertices" << std::endl;

        ASSERT_MSG(source < num_vertices, "Invalid source vertex!");
        ASSERT_MSG(target < num_vertices, "Invalid target vertex!");
      }

      if(source == target) {
        logstream(LOG_FATAL)
          << "Attempting to add self edge (" << source << " -> " << target <<  ").  "
          << "This operation is not permitted in GraphLab!" << std::endl;
        ASSERT_MSG(source != target, "Attempting to add self edge!");
      }

      // Add the edge to the set of edge data (this copies the edata)
      edge_id_type edge_id = 0;
      while (true) {
        edge_id = edges->size();
        bool ok = edges->hsetnx(Int2HexString(edge_id), "E", serialize_to_string(edata));
        if (ok) {
          break;
        }
      }
      vector<pair<string, string> > edge;
      auto source_str = Int2HexString(source);
      edge.push_back(make_pair("S", source_str));
      auto target_str = Int2HexString(target);
      edge.push_back(make_pair("T", target_str));
      auto edge_id_str = Int2HexString(edge_id);
      edges->hmset(edge_id_str, edge);

      // Add the edge id to in and out edge maps
      in_edges->hset(target_str, source_str, edge_id_str);
      out_edges->hset(source_str, target_str, edge_id_str);

      return edge_id;
    } // End of add edge

    /** \brief Returns a reference to the data stored on the vertex v. */
    void set_vertex_data(vertex_id_type v, const VertexData& data) {
      vertices->set(Int2HexString(v), serialize_to_string(data));
    } // end of data(v)

    /** \brief Returns a constant reference to the data stored on the vertex v */
    const VertexData& vertex_data(vertex_id_type v) const {
      DCHECK_LT(v, vertices->size());
      VertexData vertex;
      deserialize_from_string(vertices->get(v), vertex);
      return vertex;
    } // end of data(v)

    /** \brief Returns a reference to the data stored on the edge source->target. */
    void set_edge_data(vertex_id_type source, vertex_id_type target, const EdgeData& data) {
      DCHECK_LT(source, vertices->size());
      DCHECK_LT(target, vertices->size());
      auto ans = find(source, target);
      // We must find the edge!
      ASSERT_TRUE(ans.first);
      // the edge id should be valid!
      DCHECK_LT(ans.second, edges->size());
      set_edge_data(ans.second, data);
    } // end of edge_data(u,v)

    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const EdgeData& edge_data(vertex_id_type source, vertex_id_type target) const {
      DCHECK_LT(source, vertices->size());
      DCHECK_LT(target, vertices->size());
      auto ans = find(source, target);
      // We must find the edge!
      ASSERT_TRUE(ans.first);
      // the edge id should be valid!
      DCHECK_LT(ans.second, edges->size());
      return edge_data(ans.second);
    } // end of edge_data(u,v)

    /** \brief Returns a reference to the data stored on the edge e */
    void set_edge_data(edge_id_type edge_id, const EdgeData& data) {
      edges->hset(Int2HexString(edge_id), "E", serialize_to_string(data));
    }

    /** \brief Returns a constant reference to the data stored on the edge e */
    const EdgeData& edge_data(edge_id_type edge_id) const {
      DCHECK_LT(edge_id, edges->size());
      EdgeData e;
      deserialize_from_string(edges->hget(Int2HexString(edge_id), "E"), e);
      return e;
    }

    /** \brief Returns the source vertex of an edge. */
    vertex_id_type source(edge_id_type edge_id) const {
      //      DCHECK_LT(edge_id, edges.size());
      return HexString2Int<vertex_id_type>(edges->hget(Int2HexString(edge_id), "S"));
    }

    /** \brief Returns the destination vertex of an edge. */
    vertex_id_type target(edge_id_type edge_id) const {
      //      DCHECK_LT(edge_id, edges.size());
      return HexString2Int<vertex_id_type>(edges->hget(Int2HexString(edge_id), "T"));
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    const vertex_color_type& color(vertex_id_type vertex) const {
      DCHECK_LT(vertex, vcolors->size());
      return HexString2Int<vertex_color_type>(vcolors->get(Int2HexString(vertex)));
    }

    void set_color(vertex_id_type vid, const vertex_color_type& col) {
      vcolors->set(Int2HexString(vid), Int2HexString(col));
    }

    /** \brief This function constructs a heuristic coloring for the
        graph and returns the number of colors */
    size_t compute_coloring() {
      size_t num_vertices = num_vertices();
      // Reset the colors
      for(auto v = 0; v < num_vertices; ++v) set_color(v, 0);
      // construct a permuation of the vertices to use in the greedy
      // coloring. \todo Should probably sort by degree instead when
      // constructing greedy coloring.
      std::vector<std::pair<vertex_id_type, vertex_id_type> > permutation(num_vertices);

      for(auto v = 0; v < num_vertices; ++v)
        permutation[v] = std::make_pair(-num_in_neighbors(v), v);
      //      std::random_shuffle(permutation.begin(), permutation.end());
      std::sort(permutation.begin(), permutation.end());
      // Recolor
      size_t max_color = 0;
      for(auto i = 0; i < permutation.size(); ++i) {
        std::set<vertex_color_type> neighbor_colors;
        const auto& vid = permutation[i].second;
        // Get the neighbor colors
        for(const auto& eid : in_edge_ids(vid)){
          const auto& neighbor_vid = source(eid);
          const auto& neighbor_color = color(neighbor_vid);
          neighbor_colors.insert(neighbor_color);
        }
        for(auto const& eid : out_edge_ids(vid)){
          const auto& neighbor_vid = target(eid);
          const auto& neighbor_color = color(neighbor_vid);
          neighbor_colors.insert(neighbor_color);
        }

        vertex_color_type vertex_color = 0;
        for(const auto& neighbor_color : neighbor_colors) {
          if(vertex_color != neighbor_color) break;
          else vertex_color++;
          // Ensure no wrap around
          ASSERT_NE(vertex_color, 0);
        }
        set_color(v, vertex_color);
        max_color = std::max(max_color, size_t(vertex_color) );

      }
      // Return the NUMBER of colors
      return max_color + 1;
    } // end of compute coloring


    /**
     * \brief Check that the colors satisfy a valid coloring of the graph.
     * return true is coloring is valid;
     */
    bool valid_coloring() const {
      auto num_vertices = num_vertices();
      for(auto vid = 0; vid < num_vertices; ++vid) {
        const auto& vertex_color = color(vid);
        const auto& in_edges = in_edge_ids(vid);
        // Get the neighbor colors
        for(const auto& eid : in_edges){
          const auto& neighbor_vid = source(eid);
          const auto& neighbor_color = color(neighbor_vid);
          if(vertex_color == neighbor_color) return false;
        }
      }
      return true;
    }


    /** \brief Return the edge ids of the edges arriving at v */
    edge_list_type in_edge_ids(vertex_id_type v) const {
      DCHECK_LT(v, in_edges->size());
      edge_list_type edge_list;
      Redis::Reply reply = in_edges->hkeys(Int2HexString(v));
      for (auto i = 0; i < reply.size(); ++i) {
        edge_list.push_back(HexString2Int<edge_id_type>(reply[i]));
      }
      return edge_list;
    } // end of in edges

    /** \brief Return the edge ids of the edges leaving at v */
    edge_list_type out_edge_ids(vertex_id_type v) const {
      DCHECK_LT(v, out_edges->size());
      edge_list_type edge_list;
      Redis::Reply reply = out_edges->hkeys(Int2HexString(v));
      for (auto i = 0; i < reply.size(); ++i) {
        edge_list.push_back(HexString2Int<edge_id_type>(reply[i]));
      }
      return edge_list;
    } // end of out edges

    /** \brief Get the set of in vertices of vertex v */
    std::vector<vertex_id_type> in_vertices(vertex_id_type v) const {
      DCHECK_LT(v, in_edges->size());
      std::vector<vertex_id_type> results;
      Redis::Reply reply = in_edges->hvals(Int2HexString(v));
      for (auto i = 0; i < reply.size(); ++i) {
        results.push_back(HexString2Int<vertex_id_type>(reply[i]));
      }
      return results;
    }

    /** \brief Get the set of out vertices of vertex v */
    std::vector<vertex_id_type> out_vertices(vertex_id_type v) const {
      DCHECK_LT(v, out_edges->size());
      std::vector<vertex_id_type> results;
      Redis::Reply reply = out_edges->hvals(Int2HexString(v));
      for (auto i = 0; i < reply.size(); ++i) {
        results.push_back(HexString2Int<vertex_id_type>(reply[i]));
      }
      return results;
    }


    /** \brief count the number of times the graph was cleared and rebuilt */
    size_t get_changeid() const {
      return changeid;
    }

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
      clear();
      {
        std::vector<VertexData> vertices_;
        arc >> vertices_;
        for (int i = 0; i < vertices_.size(); ++i) {
          ASSERT_EQ(i, add_vertex(vertices_[i]));
        }
      }
      {
        std::vector<edge> edges_;
        arc >> edges_;
        for (int i = 0; i < edges_.size(); ++i) {
          ASSERT_EQ(i, add_edge(edges_[i].source(),
                                  edges_[i].target(),
                                  edges_[i].data()));
        }
      }
      {
        std::vector<std::vector<edge_id_type> > ignored;
        arc >> ignored; // in_edges
        arc >> ignored; // out_edges
      }
      {
        std::vector<vertex_color_type> colors_;
        arc >> colors_;
        for (int i = 0; i < colors_.size(); ++i) {
          set_color(i, colors_[i]);
        }
      }
      {
        bool finalized;
        arc >> finalized;
      }
    } // end of load

    /** \brief Save the graph to an archive */
    void save(oarchive& arc) const {
      assert(false);
    } // end of save


    /** \brief Load the graph from a file */
    void load(const std::string& filename) {
      std::ifstream fin(filename.c_str());
      iarchive iarc(fin);
      iarc >> *this;
      fin.close();
    } // end of load


    /**
     * \brief save the graph to the file given by the filename
     */
    void save(const std::string& filename) const {
      std::ofstream fout(filename.c_str());
      oarchive oarc(fout);
      oarc << *this;
      fout.close();
    } // end of save

    /**
     * \brief save the adjacency structure to a text file.
     *
     * Save the adjacency structure as a text file in:
     *    src_Id, dest_Id \n
     *    src_Id, dest_Id \n
     * format.
     */
    void save_adjacency(const std::string& filename) const {
      std::ofstream fout(filename.c_str());
      ASSERT_TRUE(fout.good());
      vector<string> fields;
      fields.push_back("S");
      fields.push_back("T");
      for(size_t i = 0; i < edges->size(); ++i) {
        Reply reply = edges->hmget(Int2HexString(i), fields);
        fout << HexString2Int<vertex_id_type>(reply[0]) << ", " << HexString2Int<vertex_id_type>(reply[1]) << "\n";
        ASSERT_TRUE(fout.good());
      }
      fout.close();
    }

    /**
     * builds a topological_sort of the graph returning it in topsort.
     *
     * \param[out] topsort Resultant topological sort of the graph vertices.
     *
     * function will return false if graph is not acyclic.
     */
    bool topological_sort(std::vector<vertex_id_type>& topsort) const {
      topsort.clear();
      num_vertices = num_vertices();
      topsort.reserve(num_vertices);

      std::vector<size_t> indeg;
      indeg.resize(num_vertices);
      std::queue<vertex_id_type> q;
      for (size_t i = 0;i < num_vertices; ++i) {
        indeg[i] = in_edge_ids(i).size();
        if (indeg[i] == 0) {
          q.push(i);
        }
      }

      while (!q.empty()) {
        vertex_id_type v = q.front();
        q.pop();
        topsort.push_back(v);
        for(const auto& eid : out_edge_ids(v)) {
          vertex_id_type destv = target(eid);
          --indeg[destv];
          if (indeg[destv] == 0) {
            q.push(destv);
          }
        }
      }

      if (q.empty() && topsort.size() != num_vertices) return false;
      return true;
    } // end of topological sort


  private:
    /** Internal edge class  */
    class edge {
      vertex_id_type _source;
      vertex_id_type _target;
      EdgeData _data;
    public:
      edge() : _source(-1), _target(-1) { }
      edge(const edge& other) :
        _source(other.source()), _target(other.target()),
        _data(other.data()) { }

      inline vertex_id_type source() const { return _source; }
      inline vertex_id_type target() const { return _target; }
      inline const EdgeData& data() const { return _data; }

      void load(iarchive& arc) {
        arc >> _source
            >> _target
            >> _data;
      }

      void save(oarchive& arc) const {
        arc << _source
            << _target
            << _data;
      }
    }; // end of edge_data

    // PRIVATE DATA MEMBERS ===================================================>
    /** The vertex data is simply node id -> VALUE */
    Redis *vertices;

    /** The edge data is edge id -> HASH(Source, Destination, Edge) */
    Redis *edges;

    /** The in edge data is dest id -> HASH(Source, EdgeID) */
    Redis *in_edges;

    /** The out edge data is src id -> HASH(Destination, EdgeID) */
    Redis *out_edges;

    /** The vertex colors specified by the user. **/
    /** The vcolor data is simply node id -> COLOR */
    Redis *vcolors;

    /** increments whenever the graph is cleared. Used to track the
     *  changes to the graph structure  */
    size_t changeid;
  }; // End of graph

  template <typename T>
  T HexString2Int(const string& data) {
    return coerce::as<T>(data, coerce::tag::hex());
  }

  template <typename T>
  string Int2HexString(const T& data) {
    return coerce::as<string>(data, coerce::tag::hex());
  }

  template<typename VertexData, typename EdgeData>
  std::ostream& operator<<(std::ostream& out,
                           const graph<VertexData, EdgeData>& graph) {
    // out << "Printing vertex data:\n";
    // for(size_t i = 0; i < graph.num_vertices(); ++i) {
    //   out << "V_" << i << ": D[" << graph.vertex_data(i) << "]\n";
    // }
    // out << "Printing edge data:\n";
    // for(size_t i = 0; i < graph.num_edges(); ++i) {
    //   out << "(V_" << graph.source(i) << "-> V_" << graph.target(i) << "): "
    //       << "D[" << graph.edge_data(i) << "]\n";
    // }
    // return out;
    // Print adjacency List
    typedef typename graphlab::graph<VertexData, EdgeData>::vertex_id_type
      vertex_id_type;
    typedef typename graphlab::graph<VertexData, EdgeData>::edge_id_type
      edge_id_type;
    auto num_vertices = graph.num_vertices();
    for(vertex_id_type vid = 0; vid < num_vertices; ++vid) {
      for(const auto& eid : graph.out_edge_ids(vid))
        out << vid << ", " << graph.target(eid) << '\n';
    }
    return out;
  }
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
