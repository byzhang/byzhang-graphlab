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
 * Also contains code that is Copyright 2011 Yahoo! Inc.  All rights
 * reserved.  
 *
 * Contributed under the iCLA for:
 *    Joseph Gonzalez (jegonzal@yahoo-inc.com) 
 *
 */



#ifndef GRAPHLAB_VERTEX_FUNCTOR_SET_HPP
#define GRAPHLAB_VERTEX_FUNCTOR_SET_HPP


#include <vector>

#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/parallel/pthread_tools.hpp>



namespace graphlab {

  template<typename UpdateFunctor>
  class vertex_functor_set {
  public:
    typedef UpdateFunctor update_functor_type;


  private:
    
    class vfun_type {
    private:
      update_functor_type functor;
      spinlock lock;
      bool is_set;

    public:
      vfun_type() : is_set(false) { }
      /** returns true if set for the first time */
      inline bool set(const update_functor_type& other) {
        lock.lock();
        const bool already_set(is_set);
        if(is_set) functor += other;
        else { functor = other; is_set = true; }
        lock.unlock();
        return !already_set;
      }

      /** returns true if set for the first time */
      inline bool set(const update_functor_type& other, 
                      double& ret_priority) {
        lock.lock();
        const bool already_set(is_set);
        if(is_set) functor += other;
        else { functor = other; is_set = true; }
        ret_priority = functor.priority();
        lock.unlock();
        return !already_set;
      }

      inline update_functor_type get() {
        update_functor_type ret;
        lock.lock();
        ASSERT_TRUE(is_set);
        ret = functor;
        is_set = false;
        lock.unlock();
        return functor;
      }
      
      void reset() {
        lock.lock();
        is_set = false;
        lock.unlock();
      }
      
      void reset_unsync() {
        is_set = false;
      }
      
      bool priority(double& ret_priority) const {        
        lock.lock();
        const bool was_set = is_set;
        double& priority = functor.priority();
        lock.unlock();
        return std::make_pair(was_set, priority);               
      }
      inline bool test_and_get(update_functor_type& ret) {
        lock.lock();
        const bool success(is_set);
        if(success) {
          ret = functor;
          is_set = false;
        }
        lock.unlock();
        return success;
      }
      
      inline bool read_value(update_functor_type& ret) {
        lock.lock();
        const bool success(is_set);
        if(success) {
          ret = functor;
        }
        lock.unlock();
        return success;
      }

    }; // end of vfun_type;

   
    typedef std::vector< vfun_type > vfun_set_type; 
    vfun_set_type vfun_set;


  public:
    /** Initialize the per vertex task set */
    vertex_functor_set(size_t num_vertices = 0) :
      vfun_set(num_vertices) { }

    /**
     * Resize the internal locks for a different graph
     */
    void resize(size_t num_vertices) {
      vfun_set.resize(num_vertices);
    }

    bool priority(vertex_id_type vid, double& ret_priority) const {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].priority(ret_priority);
    } // end of priority

    
    /** Add a task to the set returning false if the task was already
        present. Promote task to max(old priority, new priority) */
    bool add(const vertex_id_type& vid, 
             const update_functor_type& fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].set(fun);
    } // end of add task to set 

    
    /** Add a task to the set returning false if the task was already
        present. Promote task to max(old priority, new priority) */
    bool add(const vertex_id_type& vid, 
             const update_functor_type& fun,
             double& ret_priority) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].set(fun, ret_priority);
    } // end of add task to set 



    bool test_and_get(const vertex_id_type& vid,
                      update_functor_type& ret_fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].test_and_get(ret_fun);
    }

    bool read_value(const vertex_id_type& vid,
                      update_functor_type& ret_fun) {
      ASSERT_LT(vid, vfun_set.size());
      return vfun_set[vid].read_value(ret_fun);
    }

    size_t size() const { 
      return vfun_set.size(); 
    }
    
    void clear_unsync() {
      for (size_t i = 0; i < vfun_set.size(); ++i) vfun_set[i].reset_unsync();
    }
    
  }; // end of vertex functor set

}; // end of namespace graphlab


#endif

