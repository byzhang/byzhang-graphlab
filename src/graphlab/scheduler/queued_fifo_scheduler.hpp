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


#ifndef GRAPHLAB_QUEUED_FIFO_SCHEDULER_HPP
#define GRAPHLAB_QUEUED_FIFO_SCHEDULER_HPP

#include <algorithm>
#include <queue>


#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>


#include <graphlab/scheduler/ischeduler.hpp>
#include <graphlab/scheduler/terminator/iterminator.hpp>
#include <graphlab/scheduler/vertex_functor_set.hpp>

#include <graphlab/scheduler/terminator/critical_termination.hpp>
#include <graphlab/options/options_map.hpp>


#include <graphlab/macros_def.hpp>
namespace graphlab {

  /**
   * \ingroup group_schedulers 
   *
   * This class defines a multiple queue approximate fifo scheduler.
   * Each processor has its own in_queue which it puts new tasks in
   * and out_queue which it pulls tasks from.  Once a processors
   * in_queue gets too large, the entire queue is placed at the end of
   * the shared master queue.  Once a processors out queue is empty it
   * grabs the next out_queue from the master.
   */
  template<typename Graph, typename UpdateFunctor>
  class queued_fifo_scheduler : public ischeduler<Graph, UpdateFunctor> {
  
  public:

    typedef ischeduler<Graph, UpdateFunctor> base;
    typedef typename base::graph_type graph_type;
    typedef typename base::vertex_id_type vertex_id_type;
    typedef typename base::update_functor_type update_functor_type;

    typedef std::deque<vertex_id_type> queue_type;

  private:
    vertex_functor_set<update_functor_type> vfun_set;
    std::deque<queue_type> master_queue;
    mutex master_lock;
    size_t sub_queue_size;
    std::vector<queue_type> in_queues;
    std::vector<mutex> in_queue_locks;
    std::vector<queue_type> out_queues;
    // Terminator
    critical_termination term;

  public:

    queued_fifo_scheduler(const graph_type& graph, 
                          size_t ncpus,
                          const options_map& opts) :
      vfun_set(graph.num_vertices()), 
      sub_queue_size(100), 
      in_queues(ncpus), in_queue_locks(ncpus), 
      out_queues(ncpus), term(ncpus) { 
      opts.get_option("queuesize", sub_queue_size);
    }

    void start() { 
      master_lock.lock();
      for (size_t i = 0;i < in_queues.size(); ++i) {
        master_queue.push_back(in_queues[i]);
        in_queues[i].clear();
      }
      master_lock.unlock();
      term.reset(); 
    }

    void schedule(const vertex_id_type vid, 
                  const update_functor_type& fun) {      
      if (vfun_set.add(vid, fun)) {
        const size_t cpuid = random::rand() % in_queues.size();
        in_queue_locks[cpuid].lock();
        queue_type& queue = in_queues[cpuid];
        queue.push_back(vid);
        if(queue.size() > sub_queue_size) {
          master_lock.lock();
          queue_type emptyq;
          master_queue.push_back(emptyq);
          master_queue.back().swap(queue);
          master_lock.unlock();
        }
        in_queue_locks[cpuid].unlock();
        term.new_job(cpuid);
      } 
    } // end of schedule

    void schedule_from_execution_thread(const size_t cpuid,
                                        const vertex_id_type vid, 
                                        const update_functor_type& fun) {      
      if (vfun_set.add(vid, fun)) {
        ASSERT_LT(cpuid, in_queues.size());
        in_queue_locks[cpuid].lock();
        queue_type& queue = in_queues[cpuid];
        queue.push_back(vid);
        if(queue.size() > sub_queue_size) {
          master_lock.lock();
          queue_type emptyq;
          master_queue.push_back(emptyq);
          master_queue.back().swap(queue);
          master_lock.unlock();
        }
        in_queue_locks[cpuid].unlock();
        term.new_job(cpuid);
      } 
    } // end of schedule

    void schedule_all(const update_functor_type& fun,
                      const std::string& order) {
      if(order == "shuffle") {
        std::vector<vertex_id_type> permutation = 
          random::permutation<vertex_id_type>(vfun_set.size());       
        foreach(vertex_id_type vid, permutation)  schedule(vid, fun);
      } else {
        for (vertex_id_type vid = 0; vid < vfun_set.size(); ++vid)
          schedule(vid, fun);      
      }
    } // end of schedule_all

    void completed(const size_t cpuid,
                   const vertex_id_type vid,
                   const update_functor_type& fun) {
      term.completed_job();
    }


    sched_status::status_enum 
    get_specific(vertex_id_type vid,
                 update_functor_type& ret_fun) {
      bool get_success = vfun_set.test_and_get(vid, ret_fun); 
      if (get_success) return sched_status::NEW_TASK;
      else return sched_status::EMPTY;
    }

    /** Get the next element in the queue */
    sched_status::status_enum get_next(const size_t cpuid,
                                       vertex_id_type& ret_vid,
                                       update_functor_type& ret_fun) {
      // if the local queue is empty try to get a queue from the master
      while(1) {
        if(out_queues[cpuid].empty()) {
          master_lock.lock();
          if(!master_queue.empty()) {
            out_queues[cpuid].swap(master_queue.front());
            master_queue.pop_front();
          }
          master_lock.unlock();
        }
        // if the local queue is still empty see if there is any local
        // work left
        in_queue_locks[cpuid].lock();
        if(out_queues[cpuid].empty() && !in_queues[cpuid].empty()) {
          out_queues[cpuid].swap(in_queues[cpuid]);
        }
        in_queue_locks[cpuid].unlock();
        // end of get next
        queue_type& queue = out_queues[cpuid];
        if(!queue.empty()) {
          ret_vid = queue.front();
          queue.pop_front();
          if(vfun_set.test_and_get(ret_vid, ret_fun)) {
            return sched_status::NEW_TASK;
          }
        } else {
          return sched_status::EMPTY;
        }
      }
    } // end of get_next_task

    iterminator& terminator() { return term; }

    size_t num_joins() const {
      return vfun_set.num_joins();
    }
    /**
     * Print a help string describing the options that this scheduler
     * accepts.
     */
    static void print_options_help(std::ostream& out) { 
      out << "\t queuesize=100: the size at which a subqueue is "
          << "placed in the master queue" << std::endl;
    }


  }; 


} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif

