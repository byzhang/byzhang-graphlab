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

#ifndef GRAPHLAB_DISTRIBUTED_GRAPH_HPP
#define GRAPHLAB_DISTRIBUTED_GRAPH_HPP

#include <omp.h>
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
#include <functional>
#include <fstream>

#include <boost/lexical_cast.hpp>


#include <graphlab/logger/logger.hpp>
#include <graphlab/logger/assertions.hpp>


#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/util/mpi_tools.hpp>

#include <graphlab/serialization/iarchive.hpp>
#include <graphlab/serialization/oarchive.hpp>
#include <graphlab/graph/graph.hpp>




#include <graphlab/macros_def.hpp>
namespace graphlab { 


  // CLASS GRAPH ==============================================================>  
  template<typename VertexData, typename EdgeData>
  class distributed_graph {
  public:

    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;

    /// The type of a vertex is a simple size_t
    typedef graphlab::vertex_id_type vertex_id_type;
    /// Type for vertex colors 
    typedef graphlab::vertex_color_type vertex_color_type;


    /// The type of the local graph used to store the graph data 
    typedef graphlab::graph<VertexData, EdgeData> local_graph_type;
    /** 
     * The type of the local vertex id.  While this is the same as the
     * vertex_id_type giving it a separate name will make method calls
     * clearer.
     */
    typedef typename local_graph_type::vertex_id_type lvid_type;

    
    /** This class represents an edge with source() and target()*/
    class edge_type {
      typedef typename local_graph_type::edge_type local_edge_type;
      local_edge_type local_edge;
      const distributed_graph* graph_ptr;
      edge_type(const distributed_graph* graph_ptr, 
                const local_edge_type local_edge) : 
        graph_ptr(graph_ptr), local_edge(local_edge) { }     
    public:
      edge_type() : graph_ptr(NULL) { }
      inline vertex_id_type source() const {
        ASSERT_FALSE(empty()); 
        // TODO: implement local to global vid conversion
      }
      inline vertex_id_type target() const { 
        ASSERT_FALSE(empty());
        // TODO: implement local to global vid conversion
      }
      bool empty() const { return graph_ptr == NULL; }
      friend class distributed_graph;
    }; // end of class edge_type.
    friend class edge_type;


    /** This class represents a lazy list of edge_type. */
    class edge_list_type {
    public:
      typedef typename edge_type::local_edge_type local_edge_type;
      typedef typename local_graph_type::edge_list_type local_edge_list_type;
      typedef typename local_edge_list_type::const_iterator 
      local_const_iterator_type;
      struct edge_functor : 
        public std::unary_function<local_edge_type, edge_type> {
        const distributed_graph* graph_ptr;
        edge_functor(const distributed_graph* graph_ptr) : 
          graph_ptr(graph_ptr) { }
        edge_type operator()(const local_edge_type& edge) const {
          return edge_type(graph_ptr, edge);
        }
      }; // end of edge_functor
      typedef boost::
      transform_iterator<edge_functor, local_const_iterator_type> iterator;

      typedef iterator const_iterator;
      typedef edge_type value_type;
    private:
      iterator begin_iter, end_iter;
    public:
      edge_list_type(const distributed_graph* graph_ptr = NULL, 
                     local_const_iterator_type begin_iter = 
                     local_const_iterator_type(), 
                     local_const_iterator_type end_iter =
                     local_const_iterator_type()) :
        begin_iter(begin_iter, edge_functor(graph_ptr)),
        end_iter(end_iter, edge_functor(graph_ptr)) { }
      size_t size() const { return end_iter - begin_iter; }
      edge_type operator[](size_t i) const {
        ASSERT_LT(i, size()); return *(begin_iter + i);
      }
      const_iterator begin() const { return begin_iter; }
      const_iterator end() const { return end_iter; }
      bool empty() const { return size() == 0; }  
    }; // end of class edge_list_type


  private:
    /**
     * The vertex record stores information associated with each
     * vertex on this proc
     */
    struct vertex_record {
      /// The official owning processor for this vertex
      procid_t owner; 
      /// The local vid of this vertex on this proc
      lvid_type lvid;
      /// The set of proc that mirror this vertex.
      std::vector<procid_t> mirrors;
      vertex_record() : owner(-1), lvid(-1) { }
    }; // vertex_record

    /// The master vertex record map
    typedef boost::unordered_map<vertex_id_type, vertex_record> 
    vid2record_type;


    /// temporary map for vertexdata
    boost::unordered_map<vertex_id_type, vertex_data_type> tmp_vdata;
    mutex tmp_vdata_lock;
      
    // PRIVATE DATA MEMBERS ===================================================>    
    
    /** The rpc interface for this class */
    mutable dc_dist_object<distributed_graph> rpc;

    /** The local graph data */
    local_graph_type local_graph;
    mutex local_graph_lock;
    
    /** The map from global vertex ids to vertex records */
    vid2record_type vid2record;
    mutex vid2record_lock;
    
    /** The map from local vertex ids back to global vertex ids */
    std::vector<vertex_id_type> lvid2vid;

    /** The global number of vertices and edges */
    size_t nverts, nedges;
   
  public:

    // CONSTRUCTORS ==========================================================>
    distributed_graph(distributed_control& dc) : rpc(dc, this) {
      rpc.barrier();
    }



    // METHODS ===============================================================>

    /**
     * Finalize is used to complete graph ingress by resolving vertex
     * ownship and completing local data structures.
     */
    void finalize();
            
    /** \brief Get the number of vertices */
    size_t num_vertices() const { return nverts; }

    /** \brief Get the number of edges */
    size_t num_edges() const { return nedges; }

    /** \brief Get the number of vertices local to this proc */
    size_t num_local_vertices() const { return local_graph.num_vertices(); }

    /** \brief Get the number of edges local to this proc */
    size_t num_local_edges() const { return local_graph.num_edges(); }

    /** 
     * \brief Creates a vertex containing the vertex data
     */
    void add_vertex(const vertex_id_type& vid, 
                    const VertexData& vdata = VertexData() );

    /**
     * \brief Creates an edge connecting vertex source to vertex target.  
     */
    void add_edge(vertex_id_type source, vertex_id_type target, 
                  const EdgeData& edata = EdgeData());

  private:
    // Helper type used to synchronize the vertex data and assignments
    struct vdata_shuffle_record {
      vertex_id_type vid;
      procid_t owner;
      std::vector<procid_t> mirrors;
      vertex_data_type vdata;
      void load(iarchive& arc) { arc >> vid >> owner >> mirrors >> vdata; }
      void save(oarchive& arc) const { arc << vid << owner << mirrors << vdata; }
    };

    // HELPER ROUTINES =======================================================>    
    procid_t edge_to_proc(vertex_id_type source, vertex_id_type target) const {
      if(source > target) std::swap(source, target);
      boost::hash< std::pair<vertex_id_type, vertex_id_type> > hash_function;
      return hash_function(std::make_pair(source, target)) % rpc.numprocs();
    }

    procid_t is_local(vertex_id_type source, vertex_id_type target) const {
      return edge_to_proc(source, target) == rpc.procid();
    }

    procid_t vertex_to_init_proc(vertex_id_type vid) const { 
      return vid % rpc.numprocs();
    }    
    
    bool is_local_init(vertex_id_type vid) const {
      return vertex_to_init_proc(vid) == rpc.procid();
    }

  }; // End of graph











  //////////////////////////////////////////////////////////////////////////////
  //// Implementation 

  template<typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::finalize()  {   
    rpc.full_barrier();
    // Check conditions on graph
    ASSERT_EQ(local_graph.num_vertices(), vid2record.size());
    // resize local vid map
    lvid2vid.resize(vid2record.size());
    // For all the vertices that this processor has seen determine the
    // "negotiator" and send that machine the negotiator.
    typedef std::vector< std::vector<vertex_id_type> >  proc2vids_type;
    proc2vids_type proc2vids(rpc.numprocs());
    {
      typedef typename vid2record_type::value_type pair_type;
      foreach(const pair_type& pair, vid2record) {
        const vertex_id_type vid = pair.first;
        const vertex_record& record = pair.second;
        proc2vids[vertex_to_init_proc(vid)].push_back(vid);
        lvid2vid[record.lvid] = vid;
      }
      // The returned local vertices are the vertices from each
      // machine for which this machine is a negotiator.
      mpi_tools::all2all(proc2vids, proc2vids);
    }

    // Construct the vid2procs map by inverting the proc2vids map
    typedef boost::unordered_map<vertex_id_type, std::vector<procid_t> > 
      vid2procs_type;
    vid2procs_type vid2procs;
    for(procid_t proc = 0; proc < rpc.numprocs(); ++proc)
      for(size_t i = 0; i < proc2vids[proc].size(); ++i) 
        vid2procs[proc2vids[proc][i]].push_back(proc);
    // clear the old local vertices
    proc2vids_type().swap(proc2vids);

    // Construct the assignments
    std::vector< std::vector<vdata_shuffle_record> > 
      vdata_shuffle(rpc.numprocs());
    std::vector<size_t> counts(rpc.numprocs());
    typedef vid2procs_type::value_type pair_type;
    foreach(const pair_type& pair, vid2procs) {
      // Find the best (least loaded) processor to assign the vertex.
      std::pair<size_t, procid_t> 
        best_asg(counts[pair.second[0]], pair.second[0]);
      foreach(procid_t proc, pair.second)
        best_asg = std::min(best_asg, std::make_pair(counts[proc], proc));
      vdata_shuffle_record record;
      record.vid = pair.first;
      record.owner = best_asg.second;
      record.mirrors = pair.second;
      record.vdata = tmp_vdata[record.vid];
      counts[record.owner]++;
      foreach(procid_t proc, pair.second)
        vdata_shuffle[proc].push_back(record);
    } // end of loop over 
    foreach(size_t c, counts) std::cout << c << '\t';
    std::cout << std::endl;

    // clear the vid2procs table
    vid2procs_type().swap(vid2procs);

    // Receive assignments from coordinators
    mpi_tools::all2all(vdata_shuffle, vdata_shuffle);
    
    // Incorporate all the vertex data assigned to each machine
    for(size_t i = 0; i < vdata_shuffle.size(); ++i) {
      foreach(vdata_shuffle_record& shuffle_record, vdata_shuffle[i]) {
        ASSERT_TRUE(vid2record.find(shuffle_record.vid) != vid2record.end());
        vertex_record& record = vid2record[shuffle_record.vid];
        record.owner = shuffle_record.owner; 
        record.mirrors.swap(shuffle_record.mirrors);
        local_graph.vertex_data(record.lvid) = shuffle_record.vdata;
      } // end of loop over vdata
    } // end of loop over sending machines
   
    // std::cout << "Save debugging information" << std::endl;
    // {
    //   const std::string fname = 
    //     "file_" + boost::lexical_cast<std::string>(rpc.procid());
    //   std::ofstream fout(fname.c_str());
    //   typedef typename vid2record_type::value_type pair_type;
    //   foreach(const pair_type& pair, vid2record) {      
    //     fout << pair.first << '\t' << pair.second.owner << '\t';
    //     std::vector<bool> bitset(rpc.numprocs(), false);
    //     foreach(const procid_t& proc, pair.second.mirrors)
    //       bitset[proc] = true;
    //     for(size_t i = 0; i < bitset.size(); ++i) {
    //       fout << (bitset[i]? '1' : '0') 
    //            << (i+1 < bitset.size()? '\t' : '\n');
    //     }
    //   }
    //   fout.close();
    // }
    

  } // End of finalize
  
  
  template<typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>::
  add_vertex(const vertex_id_type& vid,
             const VertexData& vdata) {
    // determine if the vertex is local
    if(is_local_init(vid)) {
      tmp_vdata_lock.lock();
      tmp_vdata[vid] = vdata;
      tmp_vdata_lock.unlock();
    } else {
      rpc.remote_call(vertex_to_init_proc(vid),
                      &distributed_graph::add_vertex,
                      vid, vdata);
    }
  } // End of add vertex;
  
  
  template<typename VertexData, typename EdgeData>
  void distributed_graph<VertexData, EdgeData>:: 
  add_edge(vertex_id_type source, vertex_id_type target, 
           const EdgeData& edata) {
    // determine if the edge is locally managed
    if(is_local(source,target)) {
      // get (or create) the local ids for source and target 
      vid2record_lock.lock();
      // if the record was just created (it will have an invalid
      // lvid=-1) then assign it the next lvid
      vertex_record& source_rec = vid2record[source];
      if(source_rec.lvid == vertex_id_type(-1))
        source_rec.lvid = vid2record.size() - 1;
      vertex_record& target_rec = vid2record[target];
      if(target_rec.lvid == vertex_id_type(-1))
        target_rec.lvid = vid2record.size() - 1;
      vid2record_lock.unlock();
      // Add the record to the local graph
      local_graph_lock.lock();
      const vertex_id_type max_lvid = std::max(source_rec.lvid, target_rec.lvid);
      if(max_lvid >= local_graph.num_vertices()) 
        local_graph.resize(max_lvid + 1);
      local_graph.add_edge(source_rec.lvid, target_rec.lvid, edata);
      local_graph_lock.unlock();
    } else {      
      rpc.remote_call(edge_to_proc(source, target),
                      &distributed_graph::add_edge,
                      source, target, edata);
    }
  } // End of add edge


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

