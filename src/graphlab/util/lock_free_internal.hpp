#ifndef GRAPHLAB_UTIL_LOCK_FREE_INTERNAL_HPP
#define GRAPHLAB_UTIL_LOCK_FREE_INTERNAL_HPP

#include <graphlab/util/generics/integer_selector.hpp>

namespace graphlab {
namespace lock_free_internal {

template <typename index_type>
union reference_with_counter {
  struct {
    index_type val;
    index_type counter;
  } q;
  index_type& value() {
    return q.val;
  }
  index_type& counter() {
    return q.counter;
  }
  typename u_integer_selector<sizeof(index_type) * 2>::integer_type combined;
};
  
}
}
#endif
