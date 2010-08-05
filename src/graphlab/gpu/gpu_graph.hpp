/**
 * \file gpu_graph.hpp
 *
 * This file contains the template for the GPU GraphLab graph
 * data-structure.
 *
 */

#ifndef GRAPHLAB_GPU_GRAPH_HPP
#define GRAPHLAB_GPU_GRAPH_HPP

/*
#include <cassert>
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

#include <graphlab/logger/logger.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
*/

#include <graphlab/graph/graph.hpp>

#include <graphlab/gpu/gpu_errors.hpp>
#include <graphlab/gpu/gpu_stl.hpp>

#include <graphlab/macros_def.hpp>

namespace graphlab {
  
  //  template<typename VertexData, typename EdgeData> class gpu_graph;

  // CLASS GPU_GRAPH ===========================================================
  /**
   * The GPU graph data structure which behaves like a container of edge
   * and vertex data. Intuitively the graph behaves like an
   * std::vector over VertexData with data associated with select
   * pairs of indices corresponding to edge data.
   *
   * GPU memory management:
   *  - From host thread:
   *     - Use create_gpu_graph() and clone() to allocate GPU memory and
   *       create deep copies of a graph.
   *     - Use destroy_gpu_graph() to deallocate GPU memory (before
   *       calling the destructor).
   *  - From GPU threads:
   *     - Constructors create shallow copies of the graph.
   *     - Destructors assume shallow copies of the graph.
   */
  template<typename VertexData, typename EdgeData>
  class gpu_graph {

    // PUBLIC TYPES
    //==========================================================================
  public:

    //! Vertex ID type
    typedef uint32_t vertex_id_t;

    //! Edge ID type
    typedef uint32_t edge_id_t;

    /**
     * This class defines a list of edges.
     */
    class edge_list {

    public:

      class const_iterator {
      public:
        typedef edge_id_t value_type;

        const_iterator(const pair<edge_id_t, vertex_id_t>* it)
          : it(it) { }

        const edge_id_t& operator[](size_t i) const { return it[i].first; }
        const edge_id_t& operator*() const { return it->first; }
        const edge_id_t* const operator->() const {return &(it->first); }

        const_iterator& operator++() {
          ++it;
          return *this;
        }
        const_iterator operator++(int) {
          const_iterator tmp(*this);
          ++(*this);
          return tmp;
        }
        const_iterator& operator+=(size_t i) {
          it += i;
          return *this;
        }
        //! @return  Distance (this iterator) - (other iterator).
        size_t diff(const const_iterator& other) const {
          return (it - other.it);
        }

      private:
        const pair<edge_id_t, vertex_id_t>* it;
      }; // class const_iterator

      typedef const_iterator iterator; // Should not be used
      typedef edge_id_t value_type;

    private:
      const_iterator begin_it; // Points to first element
      const_iterator end_it; // One past end   

    public:
      /** Construct an empty edge list */
      edge_list() : begin_it(NULL), end_it(NULL) { }

      /** Construct an edge list from an in-memory array */
//      edge_list(const edge_id_t* begin_it, size_t len) :
//        begin_it(begin_it),  end_it(begin_it + len) { }

      /** Construct an edge list from a gpu_vector */
      edge_list(const gpu_vector<pair<edge_id_t, vertex_id_t> > vec) :
        begin_it(vec.begin()), end_it(vec.end()) { }

      /** Get the size of the edge list */
      size_t size() const { return end_it.diff(begin_it); }

      /**
       * Get the ith element.
       * WARNING: This is not bound-checked.
       */
      edge_id_t operator[](size_t i) const {
        return begin_it[i].first;
      }

      /** Get the beginning */
      const_iterator begin() const {
        return begin_it;
      }

      /** Get the end */
      const_iterator end() const {
        return end_it;
      }

    }; // End of edge_list

    /** The type of the vertex data stored in the graph */
    typedef VertexData vertex_data_type;

    /** The type of the edge data stored in the graph */
    typedef EdgeData   edge_data_type;

    // PROTECTED METHODS
    //==========================================================================
  protected:

    friend
    gpu_graph<VertexData, EdgeData>
    create_gpu_graph(const graph<VertexData, EdgeData>& g);

    // PUBLIC METHODS
    //==========================================================================
  public:

    // CONSTRUCTORS ============================================================>

    /**
     * Build an empty gpu_graph.
     */
    gpu_graph()
      : in_edges_(NULL), out_edges_(NULL), error_code_(NULL) {  }

    // METHODS =================================================================>

    /**
     * Clear all internal data
     */
    /*
    void clear() {
      vertices.clear();
      edges.clear();
      in_edges.clear();
      out_edges.clear();
    }
    */
            
    /** Get the number of vetices */
    size_t num_vertices() const {
      return vertices_.size();
    }

    /** Get the number of edges */
    size_t num_edges() const {
      return edges_.size();
    }

    /** Get the number of in edges */
    size_t num_in_neighbors(vertex_id_t v) const {
      if (v < num_vertices()) {
        return in_edges[v].size();
      } else {
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
        return (size_t)(-1);
      }
    }
    
    /** get the number of out edges */
    size_t num_out_neighbors(vertex_id_t v) const  {
      if (v < num_vertices()) {
        return out_edges[v].size();
      } else {
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
        return (size_t)(-1);
      }
    }

    /**
     * Find an edge.
     * This does bound checking.
     */
    pair<bool, edge_id_t>
    find(vertex_id_t source, vertex_id_t target) const {
      if (source >= num_vertices() || target >= num_vertices()) {
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
        return pair<bool,edge_id_t>(false,-1);
      }
      const gpu_vector<edge_id_t> source_out_edges(out_edges_[source]);
      const gpu_vector<edge_id_t> target_in_edges(in_edges_[target]);
      // O( log degree ) search
      if(target_in_edges.size() < source_out_edges.size())
        return binary_search(target_in_edges, source);
      else
        return binary_search(source_out_edges, target);
    }

    /**
     * Get the ID for the edge.
     * NOTE: It is better to use find().
     */
    edge_id_t edge_id(vertex_id_t source, vertex_id_t target) const {
      pair<bool, edge_id_t> res(find(source, target));
      if (!res.first)
        error_code_[0] = gpu_errors::ELEMENT_NOT_FOUND;
      return res.second;
    }

    /**
     * Get the edge id of the reverse of the given edge.
     * NOTE: It is better to use find().
     */
    edge_id_t rev_edge_id(edge_id_t eid) const {
      if (eid < edges_.size()) {
        edge e(edges[eid].first);
        return edge_id(e.target(), e.source());
      } else {
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
        return (size_t)(-1);
      }
    }

    /** Get the vertex data */
    VertexData& vertex_data(vertex_id_t v) {
      if (v >= vertices_.size())
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
      return vertices_[v];
    }

    /** Get the vertex data */
    const VertexData& vertex_data(vertex_id_t v) const {
      if (v >= vertices_.size())
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
      return vertices_[v];
    }

    /** Get the edge_data */
    EdgeData& edge_data(vertex_id_t source, vertex_id_t target) {
      std::pair<bool, edge_id_t> ans(find(source, target));
      if (!ans.first)
        error_code_[0] = gpu_errors::ELEMENT_NOT_FOUND;
      return edges_[ans.second].second;
    }

    /** Get the edge_data */
    const EdgeData& edge_data(vertex_id_t source, vertex_id_t target) const {
      std::pair<bool, edge_id_t> ans(find(source, target));
      if (!ans.first)
        error_code_[0] = gpu_errors::ELEMENT_NOT_FOUND;
      return edges_[ans.second].second;
    }

    /** Get the edge_data */
    EdgeData& edge_data(edge_id_t edge_id) {
      if (edge_id >= edges_.size())
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
      return edges_[edge_id].second;
    }

    /** Get the edge_data */
    const EdgeData& edge_data(edge_id_t edge_id) const {
      if (edge_id >= edges_.size())
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
      return edges_[edge_id].second;
    }

    /** get the source of the edge */
    vertex_id_t source(edge_id_t edge_id) const {
      if (edge_id >= edges_.size())
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
      return edges_[edge_id].first.source();
    }

    /** get the dest of the edge */
    vertex_id_t target(edge_id_t edge_id) const {
      if (edge_id >= edges_.size())
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
      return edges_[edge_id].first.target();    
    }

    /** Get the ids of the in edges */
    edge_list in_edge_ids(vertex_id_t v) const {
      if (v >= vertices_.size())
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
      return edge_list(in_edges_[v]);
    }

    /** Get the ids of the out edges */
    edge_list out_edge_ids(vertex_id_t v) const {
      if (v >= vertices_.size())
        error_code_[0] = gpu_errors::OUT_OF_BOUNDS;
      return edge_list(out_edges_[v]);
    }

  protected:    
    // PROTECTED TYPES ========================================================>

    /** Internal edge class  */   
    class edge {

      vertex_id_t _source;
      vertex_id_t _target;

    public:
      edge() : _source(-1), _target(-1) { }
      edge(const edge& other) :
        _source(other.source()), _target(other.target()) { }
      edge(vertex_id_t source, vertex_id_t target) :
        _source(source), _target(target)  { }

      bool operator<(const edge& other) const {
        return (_source < other._source) || 
          (_source == other._source && _target < other._target); 
      }

      inline vertex_id_t source() const { return _source; }
      inline vertex_id_t target() const { return _target; }   

    }; // end of edge_data

    // PROTECTED DATA MEMBERS =================================================>

    /** Vector of vertex data */
    gpu_vector<VertexData> vertices_;

    //! Vector of edges < <source, destination>, edge data>.
    gpu_vector<pair<edge, EdgeData> > edges_;

    /**
     * Map: target vertex -> vector of <in edge ID, source vertex> pairs.
     * Each vector must be ordered by increasing source vertex ID.
     */
    gpu_vector<pair<edge_id_t, vertex_id_t> >* in_edges_;

    /**
     * Map: source vertex -> vector of <out edge ID, target vertex> pairs.
     * Each vector must be ordered by increasing target vertex ID.
     */
    gpu_vector<pair<edge_id_t, vertex_id_t> >* out_edges_;

    //! Error code (so the CPU can check to see if a problem occured)
    gpu_errors::gpu_error_type* error_code_;

    // PROTECTED HELPERS ======================================================>

    /**
     * Find the edge in the vector in O(log length(vec)) time.
     *
     * @return pair<indicator for whether edge was found, edge ID>
     *         If the edge is not found, returns <false, size_t(-1)>.
     */
    pair<bool, edge_id_t>
    binary_search(const gpu_vector<pair<edge_id_t, vertex_id_t> >& vec,
                  vertex_id_t v) const {
      // Compare to the middle of the list
      size_t first = 0;
      size_t last = vec.size();
      while(first < last) {
        size_t mid = (first+last)/2;
        vertex_id_t mid_v = vec[mid].second;
        if(mid_v == v) // edge found
          return make_pair(true, vec[mid].first);
        // otherwise search further
        if(v < mid_v)
          last = mid; // Search left
        else // Search right
          first = mid + 1;
      }
      // The edge was not found.
      return make_pair(false,(size_t)(-1));
    } // end of binary search 

  }; // End of gpu_graph

  /**
   * Construct a gpu_graph from a CPU graph.
   */
  template<typename VertexData, typename EdgeData>
  gpu_graph<VertexData, EdgeData>
  create_gpu_graph(const graph<VertexData, EdgeData>& g) {

    typedef typename gpu_graph<VertexData, EdgeData>::edge edge;
    typedef typename gpu_graph<VertexData, EdgeData>::edge_list edge_list;
    gpu_graph<VertexData, EdgeData> gpug;

    // Create vertices_ and edges_.
    gpug.vertices_ =
      gpu_vector<VertexData>
      (g.num_vertices(), gpu_alloc<VertexData>(g.num_vertices()));
    gpug.edges_ =
      gpu_vector<pair<edge, EdgeData> >
      (g.num_edges(), gpu_alloc<pair<edge, EdgeData> >(g.num_edges()));

    // Create in_edges_ and out_edges_.
    gpu_vector<pair<edge_id_t,vertex_id_t> >* h_inout_edges =
      new gpu_vector<pair<edge_id_t,vertex_id_t> >[g.num_vertices()];

    gpug.in_edges_ =
      gpu_alloc<gpu_vector<pair<edge_id_t,vertex_id_t> > >(g.num_vertices());
    for (size_t i(0); i < g.num_vertices(); ++i) {
      h_inout_edges[i] =
        gpu_vector<pair<edge_id_t,vertex_id_t> >
        (g.num_in_neighbors(i),
         gpu_alloc<pair<edge_id_t,vertex_id_t> >(g.num_in_neighbors(i)));
      pair<edge_id_t,vertex_id_t>* h_edge_vec =
        new pair<edge_id_t,vertex_id_t>[g.num_in_neighbors(i)];
      edge_list elist(g.in_edge_ids(i));
      for (size_t j(0); j < elist.size(); ++j)
        h_edge_vec[j] = make_pair(elist[j], g.source(elist[j]));
      copy_host_to_device(h_edge_vec, elist.size(), h_inout_edges[i].begin());
      delete [] h_edge_vec;
      h_edge_vec = NULL;
    }
    copy_host_to_device(h_inout_edges, g.num_vertices(), gpug.in_edges_);

    gpug.out_edges_ =
      gpu_alloc<gpu_vector<pair<edge_id_t,vertex_id_t> > >(g.num_vertices());
    for (size_t i(0); i < g.num_vertices(); ++i) {
      h_inout_edges[i] =
        gpu_vector<pair<edge_id_t,vertex_id_t> >
        (g.num_out_neighbors(i),
         gpu_alloc<pair<edge_id_t,vertex_id_t> >(g.num_out_neighbors(i)));
      pair<edge_id_t,vertex_id_t>* h_edge_vec =
        new pair<edge_id_t,vertex_id_t>[g.num_out_neighbors(i)];
      edge_list elist(g.out_edge_ids(i));
      for (size_t j(0); j < elist.size(); ++j)
        h_edge_vec[j] = make_pair(elist[j], g.target(elist[j]));
      copy_host_to_device(h_edge_vec, elist.size(), h_inout_edges[i].begin());
      delete [] h_edge_vec;
      h_edge_vec = NULL;
    }
    copy_host_to_device(h_inout_edges, g.num_vertices(), gpug.out_edges_);

    delete [] h_inout_edges;
    h_inout_edges = NULL;

    // Create error_code_.
    gpug.error_code_ = gpu_alloc<gpu_errors::gpu_error_type>(1);

    return gpug;
  } // create_gpu_graph()

  template<typename VertexData, typename EdgeData>
  void
  destroy_gpu_graph(const gpu_graph<VertexData, EdgeData>& g) {
    size_t nverts(g.num_vertices());

    // Free vertices_, edges_.
    vertices_.free();
    edges_.free();

    // Free in_edges_, out_edges_.
    gpu_vector<pair<edge_id_t,vertex_id_t> >* h_inout_edges =
      new gpu_vector<pair<edge_id_t,vertex_id_t> >[nverts];
    copy_device_to_host(g.in_edges_, nverts, h_inout_edges);
    for (size_t i(0); i < nverts; ++i)
      h_inout_edges[i].free();
    gpu_free(g.in_edges_);
    g.in_edges_ = NULL;

    copy_device_to_host(g.out_edges_, nverts, h_inout_edges);
    for (size_t i(0); i < nverts; ++i)
      h_inout_edges[i].free();
    gpu_free(g.out_edges_);
    g.out_edges_ = NULL;

    delete [] h_inout_edges;
    h_inout_edges = NULL;

    // Free error_code_.
    gpu_free(g.error_code_);
    g.error_code_ = NULL;
  } // destroy_gpu_graph()

} // end of namespace graphlab

#include <graphlab/macros_undef.hpp>

#endif // GRAPHLAB_GPU_GRAPH_HPP
