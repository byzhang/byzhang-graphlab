#ifndef GRAPHLAB_RESIZING_COUNTING_SINK
#define GRAPHLAB_RESIZING_COUNTING_SINK

#include <graphlab/util/charstream.hpp>

namespace graphlab {

  typedef charstream_impl::resizing_array_sink<false> resizing_array_sink;
  
  /**
  Wraps a resizing array sink.
  */
  class resizing_array_sink_ref {
   private:
    resizing_array_sink* ras;
   public:
   

    typedef resizing_array_sink::char_type char_type;
    typedef resizing_array_sink::category category;

    inline resizing_array_sink_ref(resizing_array_sink& ref): ras(&ref) { }
  
    inline resizing_array_sink_ref(const resizing_array_sink_ref& other) :
      ras(other.ras) { }

    inline size_t size() const { return ras->size(); }
    inline char* c_str() { return ras->c_str(); }

    inline void clear() { ras->clear(); }
    /** the optimal buffer size is 0. */
    inline std::streamsize optimal_buffer_size() const { 
      return ras->optimal_buffer_size(); 
    }
    
    inline std::streamsize write(const char* s, std::streamsize n) {
      return ras->write(s, n);
    }
  };
  
}
#endif
