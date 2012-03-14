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
#include <boost/coerce/tag/base.hpp>
#include <boost/unordered_set.hpp>

#include <hiredispp.h>

#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>
#include <graphlab/serialization/serialization_includes.hpp>
#include <graphlab/util/random.hpp>

#include <graphlab/macros_def.hpp>

using namespace boost;
using namespace hiredispp;

namespace graphlab {
  template <typename T>
  T HexString2Int(const std::string& data) {
    return coerce::as<T>(data, coerce::tag::hex());
  }

  template <typename T>
  std::string Int2HexString(const T& data) {
    return coerce::as<std::string>(data, coerce::tag::hex());
  }

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
     the graphs to 4 billion _vertices it also helps reduce the storage
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
    typedef graphlab::vertex_id_type vertex_id_type;

    /// The type of an edge id
    typedef graphlab::edge_id_type edge_id_type;

    /// Type for vertex colors
    typedef graphlab::vertex_color_type vertex_color_type;

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;

    /* ----------------------------------------------------------------------------- */
    /* Helper data field and structures: redis_server, edge_type                     */
    /* ----------------------------------------------------------------------------- */
    struct redis_server {
      const char *host;
      int port;
      int db;
      redis_server(const char *host_="127.0.0.1", int port_=6379, int db_=0) {
        host = host_;
        port = port_;
        db = db_;
      }
    };

    class reference_vertex_data_type : public vertex_data_type {
      public:
        void set_id(const vertex_id_type& id_) {
          id = id_;
        }

        void set_graph(redis_graph* graph_) {
          graph = graph_;
        }

        virtual ~reference_vertex_data_type() {
          graph->commit_vertex_data(id, *this);
        }
      private:
        vertex_id_type id;
        redis_graph* graph;
    };

    typedef const vertex_data_type const_reference_vertex_data_type;

    class reference_edge_data_type : public edge_data_type {
      public:
        void set_id(const edge_id_type& id_) {
          id = id_;
        }

        void set_graph(redis_graph* graph_) {
          graph = graph_;
        }

        virtual ~reference_edge_data_type() {
          graph->commit_edge_data(id, *this);
        }
      private:
        edge_id_type id;
        redis_graph* graph;
    };

    typedef const edge_data_type const_reference_edge_data_type;

    // A class of edge information. Used as value type of the edge_list.
    class edge_type {
    public:
      edge_type () : _source(-1), _target(-1), _edge_id(-1), _empty(true) { }
      edge_type (const vertex_id_type _source, const vertex_id_type _target,
                 const edge_id_type _eid) :
        _source(_source), _target(_target), _edge_id(_eid), _empty(false) { }
      edge_type (const redis_graph& graph, const vertex_id_type _source, const vertex_id_type _target) :
        _source(_source)
        , _target(_target)
        , _edge_id(-1)
        , _empty(true) {
          auto result = graph._out_edges->hget(Int2HexString(_source), Int2HexString(_target));
          if (result == Redis::Nil) {
            return;
          }
          _edge_id = HexString2Int<edge_id_type>(result);
          _empty = false;
        }
      edge_type (const redis_graph& graph, const edge_id_type _eid) :
        _source(-1)
        , _target(-1)
        , _edge_id(_eid)
        , _empty(true) {
          std::vector<std::string> fields;
          fields.push_back("S");
          fields.push_back("T");
          auto reply = graph._edges->hmget(Int2HexString(_eid), fields);
          const auto& result = reply[0];
          if (result.isNil()) {
            return;
          }
          _source = HexString2Int<vertex_id_type>(result);
          _target = HexString2Int<vertex_id_type>(reply[1]);
          _empty = false;
        }
    public:
      inline vertex_id_type source() const {
        // ASSERT_FALSE(empty());
        return _source;
      }

      inline vertex_id_type target() const {
        // ASSERT_FALSE(empty());
        return _target;
      }

      inline edge_id_type edge() const {
        return _edge_id;
      }

      inline bool empty() const { return _empty; }
      // Data fields.
    private:
      vertex_id_type _source;
      vertex_id_type _target;
      edge_id_type _edge_id;
      bool _empty;
    }; // end of class edge_type.

    /** The type of the edge list */
    typedef std::vector<edge_type> edge_list_type;
    typedef std::vector<edge_type> edge_list;

  public:
    // CONSTRUCTORS ============================================================>
    /**
     * Build a basic graph
     */
    redis_graph()
      : _vertices(NULL)
      , _vcolors(NULL)
      , _edges(NULL)
      , _in_edges(NULL)
      , _out_edges(NULL)
    {
    }

    redis_graph(const redis_graph& g)
      : _vertices(NULL)
      , _vcolors(NULL)
      , _edges(NULL)
      , _in_edges(NULL)
      , _out_edges(NULL)
      , _vertex_server(g._vertex_server)
      , _color_server(g._color_server)
      , _edge_server(g._edge_server)
      , _in_edge_server(g._in_edge_server)
      , _in_edge_server(g._out_edge_server)
    {
      init(_vertex_server, _edge_server, _in_edge_server, _out_edge_server, _color_server);
    }

    void init(const redis_server& vertex_server,
              const redis_server& edge_server,
              const redis_server& in_edge_server,
              const redis_server& out_edge_server,
              const redis_server& color_server) {
      // TODO: check/handle the same server address.
      _vertex_server = vertex_server;
      _color_server = color_server;
      _edge_server = edge_server;
      _in_edge_server = in_edge_server;
      _out_edge_server = out_edge_server;
      _vertices = new Redis(vertex_server.host, vertex_server.port);
      _vcolors = new Redis(color_server.host, color_server.port);
      _edges = new Redis(edge_server.host, edge_server.port);
      _in_edges = new Redis(in_edge_server.host, in_edge_server.port);
      _out_edges = new Redis(out_edge_server.host, out_edge_server.port);
      _vertices->select(vertex_server.db);
      _vcolors->select(color_server.db);
      _edges->select(edge_server.db);
      _in_edges->select(in_edge_server.db);
      _out_edges->select(out_edge_server.db);
    }

    virtual ~redis_graph() {
      delete _vertices;
      delete _vcolors;
      delete _edges;
      delete _in_edges;
      delete _out_edges;
    }

    // METHODS =================================================================>
    /**
     * \brief Resets the graph state.
     */
    void clear() {
      _vertices->flush();
      _vcolors->flush();
      _edges->flush();
      _in_edges->flush();
      _out_edges->flush();
    }

    void finalize() {
    }

    /** \brief Get the number of _vertices */
    size_t num_vertices() const {
      if (_vertices == NULL) {
        return 0;
      }
      return _vertices->size();
    } // end of num _vertices

    /** \brief Get the number of _vertices local to this machine */
    size_t local_vertices() const {
      return num_vertices();
    } // end of num _vertices

    /** \brief Get the number of _edges */
    size_t num_edges() const {
      return _edges->size();
    } // end of num _edges

    /** \brief Get the number of in _edges of a particular vertex */
    size_t num_in_edges(vertex_id_type v) const {
      return _in_edges->hlen(Int2HexString(v));
    } // end of num _vertices

    /** \brief Get the number of out _edges of a particular vertex */
    size_t num_out_edges(vertex_id_type v) const  {
      return _out_edges->hlen(Int2HexString(v));
    } // end of num _vertices

    /** \brief Finds an edge.
        The value of the first element of the pair will be true if an
        edge from src to target is found and false otherwise. If the
        edge is found, the edge ID is returned in the second element of the pair. */
    edge_type find(vertex_id_type source, vertex_id_type target) const {
      return edge_type(*this, source, target);
    } // end of find

    edge_type reverse_edge(const edge_type& edge) const {
      return find(edge.target(), edge.source());
    }

    /**
     * \brief Creates a vertex containing the vertex data and returns the id
     * of the new vertex id. Vertex ids are assigned in increasing order with
     * the first vertex having id 0.
     */
    vertex_id_type add_vertex(const VertexData& vdata = VertexData() ) {
      vertex_id_type id = 0;
      while (true) {
        id = num_vertices();
        bool ok = _vertices->setnx(Int2HexString(id), serialize_to_string(vdata));
        if (ok) {
          break;
        }
      }
      ASSERT_EQ(1, _vcolors->setnx(Int2HexString(id), serialize_to_string(0)));
      return id;
    } // End of add vertex;

    void resize(size_t num_vertices) {
      auto current_num_vertices = this->num_vertices();
      for (size_t i = 0; i < num_vertices - current_num_vertices; ++i) {
        add_vertex();
      }
    }

    /**
     * \brief Creates an edge connecting vertex source to vertex target.  Any
     * existing data will be cleared.
     */
    edge_id_type add_edge(vertex_id_type source, vertex_id_type target,
                          const EdgeData& edata = EdgeData()) {
      auto num_vertices = this->num_vertices();
      if (source >= num_vertices || target >= num_vertices) {
        logstream(LOG_FATAL)
          << "Attempting add_edge (" << source
          << " -> " << target
          << ") when there are only " << _vertices->size()
          << " _vertices" << std::endl;

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
        edge_id = num_edges();
        bool ok = _edges->hsetnx(Int2HexString(edge_id), "E", serialize_to_string(edata));
        if (ok) {
          break;
        }
      }
      std::vector<std::pair<std::string, std::string> > edge;
      auto source_str = Int2HexString(source);
      edge.push_back(make_pair("S", source_str));
      auto target_str = Int2HexString(target);
      edge.push_back(make_pair("T", target_str));
      auto edge_id_str = Int2HexString(edge_id);
      _edges->hmset(edge_id_str, edge);
      // Add the edge id to in and out edge maps
      _in_edges->hset(target_str, source_str, edge_id_str);
      _out_edges->hset(source_str, target_str, edge_id_str);

      return edge_id;
    } // End of add edge

    /** \brief Returns a reference to the data stored on the vertex v. */
    reference_vertex_data_type vertex_data(vertex_id_type v) {
      DCHECK_LT(v, _vertices->size());
      reference_vertex_data_type vertex;
      vertex_data_type* real_vertex = &vertex;
      deserialize_from_string(_vertices->get(Int2HexString(v)), *real_vertex);
      vertex.set_id(v);
      vertex.set_graph(this);
      return vertex;
    } // end of data(v)
    void commit_vertex_data(vertex_id_type v, const VertexData& data) {
      _vertices->set(Int2HexString(v), serialize_to_string(data));
    } // end of data(v)

    /** \brief Returns a constant reference to the data stored on the vertex v */
    const_reference_vertex_data_type vertex_data(vertex_id_type v) const {
      DCHECK_LT(v, _vertices->size());
      VertexData vertex;
      deserialize_from_string(_vertices->get(Int2HexString(v)), vertex);
      return vertex;
    } // end of data(v)

    /** \brief Returns a reference to the data stored on the edge source->target. */
    reference_edge_data_type edge_data(vertex_id_type source, vertex_id_type target) {
      DCHECK_LT(source, _vertices->size());
      DCHECK_LT(target, _vertices->size());
      return edge_data(find(source, target));
    }
    reference_edge_data_type edge_data(const edge_type& edge) {
      ASSERT_FALSE(edge.empty());
      // the edge id should be valid!
      DCHECK_LT(edge.edge(), _edges->size());
      reference_edge_data_type e;
      edge_data_type* real_edge = &e;
      deserialize_from_string(_edges->hget(Int2HexString(edge.edge()), "E"), *real_edge);
      e.set_id(edge.edge());
      e.set_graph(this);
      return e;
    }
    void commit_edge_data(edge_id_type edge_id, const EdgeData& data) {
      _edges->hset(Int2HexString(edge_id), "E", serialize_to_string(data));
    }

    /** \brief Returns a constant reference to the data stored on the
        edge source->target */
    const_reference_edge_data_type edge_data(vertex_id_type source, vertex_id_type target) const {
      DCHECK_LT(source, _vertices->size());
      DCHECK_LT(target, _vertices->size());
      return edge_data(find(source, target));
    } // end of edge_data(u,v)

    /** \brief Returns a constant reference to the data stored on the edge e */
    const_reference_edge_data_type edge_data(const edge_type& edge) const {
      ASSERT_FALSE(edge.empty());
      // the edge id should be valid!
      DCHECK_LT(edge.edge(), _edges->size());
      EdgeData e;
      deserialize_from_string(_edges->hget(Int2HexString(edge.edge()), "E"), e);
      return e;
    }

    /** \brief Returns the source vertex of an edge. */
    vertex_id_type source(edge_id_type edge_id) const {
      //      DCHECK_LT(edge_id, _edges.size());
      return HexString2Int<vertex_id_type>(_edges->hget(Int2HexString(edge_id), "S"));
    }

    /** \brief Returns the destination vertex of an edge. */
    vertex_id_type target(edge_id_type edge_id) const {
      //      DCHECK_LT(edge_id, _edges.size());
      return HexString2Int<vertex_id_type>(_edges->hget(Int2HexString(edge_id), "T"));
    }

    /** \brief Returns the vertex color of a vertex.
        Only valid if compute_coloring() is called first.*/
    vertex_color_type color(vertex_id_type vertex) const {
      DCHECK_LT(vertex, _vcolors->size());
      return HexString2Int<vertex_color_type>(_vcolors->get(Int2HexString(vertex)));
    }

    void set_color(vertex_id_type vid, const vertex_color_type& col) {
      _vcolors->set(Int2HexString(vid), Int2HexString(col));
    }

    /** \brief Return the edge ids of the _edges arriving at v */
    edge_list_type in_edges(vertex_id_type v) const {
      edge_list_type edge_list;
      Redis::Reply reply = _in_edges->hgetall(Int2HexString(v));
      for (size_t i = 0; i < reply.size(); i+=2) {
        edge_list.push_back(edge_type(HexString2Int<edge_id_type>(reply[i]), v, HexString2Int<edge_id_type>(reply[i+1])));
      }
      return edge_list;
    } // end of in _edges

    /** \brief Return the edge ids of the _edges leaving at v */
    edge_list_type out_edges(vertex_id_type v) const {
      edge_list_type edge_list;
      Redis::Reply reply = _out_edges->hgetall(Int2HexString(v));
      for (size_t i = 0; i < reply.size(); i+=2) {
        edge_list.push_back(edge_type(v, HexString2Int<edge_id_type>(reply[i]), HexString2Int<edge_id_type>(reply[i+1])));
      }
      return edge_list;
    } // end of out _edges

    /** \brief Load the graph from an archive */
    void load(iarchive& arc) {
      clear();
      {
        std::vector<VertexData> vertices_;
        arc >> vertices_;
        for (size_t i = 0; i < vertices_.size(); ++i) {
          ASSERT_EQ(i, add_vertex(vertices_[i]));
        }
      }
      {
        std::vector<edge> edges_;
        arc >> edges_;
        for (size_t i = 0; i < edges_.size(); ++i) {
          ASSERT_EQ(i, add_edge(edges_[i].source(),
                                edges_[i].target(),
                                edges_[i].data()));
        }
      }
      {
        std::vector<vertex_color_type> colors_;
        arc >> colors_;
        for (size_t i = 0; i < colors_.size(); ++i) {
          set_color(i, colors_[i]);
        }
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
      std::vector<std::string> fields;
      fields.push_back("S");
      fields.push_back("T");
      for(size_t i = 0; i < _edges->size(); ++i) {
        Redis::Reply reply = _edges->hmget(Int2HexString(i), fields);
        fout << HexString2Int<vertex_id_type>(reply[0]) << ", " << HexString2Int<vertex_id_type>(reply[1]) << "\n";
        ASSERT_TRUE(fout.good());
      }
      fout.close();
    }


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
    }; // end of edge

    // PRIVATE DATA MEMBERS ===================================================>
    /** The vertex data is simply node id -> VALUE */
    Redis *_vertices;

    /** The vertex colors specified by the user. **/
    /** The vcolor data is simply node id -> COLOR */
    Redis *_vcolors;

    /** The edge data is edge id -> HASH(Source, Destination, Edge) */
    Redis *_edges;

    /** The in edge data is dest id -> HASH(Source, EdgeID) */
    Redis *_in_edges;

    /** The out edge data is src id -> HASH(Destination, EdgeID) */
    Redis *_out_edges;

    redis_server _vertex_server, _color_server, _edge_server, _in_edge_server, _out_edge_server;
  }; // End of graph

} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
