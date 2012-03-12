#ifndef GRAPHLAB_PARALLEL_LOCKFREE_PUSHBACK_HPP
#define GRAPHLAB_PARALLEL_LOCKFREE_PUSHBACK_HPP
#include <graphlab/parallel/atomic.hpp>

namespace graphlab {

namespace lockfree_push_back_impl {
  struct idx_ref {
    idx_ref(): reference_count(0), idx(0) { }
    idx_ref(size_t idx): reference_count(0), idx(idx) { }
    
    volatile int reference_count;
    atomic<size_t> idx;
    enum {
      MAX_REF = 65536
    };
    
    inline void inc_ref() {
      while (1) {
        int curref = reference_count;
        if ((curref & MAX_REF) == 0 &&
            atomic_compare_and_swap(reference_count, curref, curref + 1)) {
          break;
        }
      }      
    }

    inline void wait_till_no_ref() {
      while((reference_count & (MAX_REF - 1)) != 0);
    }
    
    inline void dec_ref() {
      __sync_fetch_and_sub(&reference_count, 1);
    }

    inline void flag_ref() {
      __sync_fetch_and_xor(&reference_count, MAX_REF);
    }
    
    inline size_t inc_idx() {
      return idx.inc_ret_last();
    }
  };
} // lockfree_push_back_impl
  
/**
 * Provides a lock free way to insert elements to the end
 * of a container. Container must provide 3 functions.
 *  - T& operator[](size_t idx)
 *  - void resize(size_t len)
 *  - size_t size()
 *
 * resize(n) must guarantee that size() >= n.
 * T& operator[](size_t idx) must succeed for idx < size() and must be 
 * safely executeable in parallel.
 * size() must be safely executeable in parallel with resize().
 */
template <typename Container, typename T = typename Container::value_type>
class lockfree_push_back {
  private:
    Container& container;
    lockfree_push_back_impl::idx_ref cur;
    mutex mut;
    float scalefactor;
  public:
    lockfree_push_back(Container& container, size_t startidx, float scalefactor = 2):
                            container(container),cur(startidx), scalefactor(scalefactor) { }

    size_t size() const {
      return cur.idx.value;
    }
    size_t push_back(T& t) {
      size_t putpos = cur.inc_idx();
      while(1) {
        cur.inc_ref();
        if (putpos < container.size()) {
          container[putpos] = t;
          cur.dec_ref();
          break;
        }
        else {
          cur.dec_ref();

          if (mut.try_lock()) {
            // ok. we need to resize
            // flag the reference and wait till there are no more references
            cur.flag_ref();
            cur.wait_till_no_ref();
            // we are exclusive here. resize
            if (putpos >= container.size()) {
              container.resize(std::max<size_t>(putpos + 1, container.size() * scalefactor));
            }
            container[putpos] = t;
            cur.flag_ref();
            mut.unlock();
            break;
          }
        }
      }
      return putpos;
    }
};

} // namespace graphlab
#endif