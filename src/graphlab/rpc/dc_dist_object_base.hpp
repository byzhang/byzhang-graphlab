#ifndef GRAPHLAB_DC_DIST_OBJECT_BASE_HPP
#define GRAPHLAB_DC_DIST_OBJECT_BASE_HPP
#include <vector>
#include <graphlab/rpc/dc_internal_types.hpp>
namespace graphlab {

namespace dc_impl {
/**
Provides some functions which allow the 
the rmi object to be modified. These functions
however, should not be called directly and are meant for internal
use. 
\todo These really should be friended, but it gets complicated due to the 
      large number of template friends
*/
class dc_dist_object_base{
 public:
  virtual void inc_calls_sent(procid_t source) = 0;
  virtual void inc_calls_received(procid_t dest) = 0;
  virtual void inc_bytes_sent(procid_t target, size_t bytes) = 0;

  virtual size_t calls_received() const = 0;
  virtual size_t calls_sent() const = 0;
};

}
}

#endif