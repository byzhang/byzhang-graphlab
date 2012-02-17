#ifndef EVENT_LOG_HPP
#define EVENT_LOG_HPP
#include <iostream>
#include <boost/bind.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/dense_bitset.hpp>
namespace graphlab {

#define EVENT_MAX_COUNTERS 256

class event_log{
 public:
  enum event_print_type{
    NUMBER,
    DESCRIPTION,
    RATE_BAR
  };
 private:
  std::ostream* out;
  volatile size_t flush_interval;
  event_print_type print_method;
  volatile bool finished;
  
  mutex m;
  conditional cond;
  
  thread printing_thread;
  std::string descriptions[EVENT_MAX_COUNTERS];
  size_t maxcounter[EVENT_MAX_COUNTERS];
  atomic<size_t> counters[EVENT_MAX_COUNTERS];
  
  timer ti;
  double prevtime;
  bool hasevents;
  size_t noeventctr;

  size_t max_desc_length;
  fixed_dense_bitset<EVENT_MAX_COUNTERS> hascounter;
 public:
  inline event_log():out(NULL),
                     flush_interval(0),
                     print_method(DESCRIPTION),
                     finished(false),
                     prevtime(0),
                     hasevents(false),
                     noeventctr(0),
                     max_desc_length(0) {
    hascounter.clear();
    for (size_t i = 0;i < EVENT_MAX_COUNTERS; ++i) {
      maxcounter[i] = 0;
    }
    printing_thread.launch(boost::bind(&event_log::thread_loop, this));
  }
  
  void initialize(std::ostream &ostrm,
                  size_t flush_interval_ms,
                  event_print_type event_print);

  void close();
  
  void thread_loop();

  void add_event_type(unsigned char eventid, std::string description);

  inline void accumulate_event(unsigned char eventid,
                               size_t count)  __attribute__((always_inline)) {
    hasevents = true;
    counters[eventid].inc(count);
  }

  void flush();

  ~event_log();
};

}
/**
 * DECLARE_EVENT_LOG(name)
 * creates an event log with a given name. This creates a variable
 * called "name" which is of type event_log. and is equivalent to:
 *
 * graphlab::event_log name;
 *
 * The primary reason to use this macro instead of just writing
 * the code above directly, is that the macro is ignored and compiles
 * to nothing when event logs are disabled.
 *
 *
 * INITIALIZE_EVENT_LOG(name, ostrm, flush_interval, printdesc)
 * ostrm is the output std::ostream object. A pointer of the stream
 * is taken so the stream should not be destroyed until the event log is closed
 * flush_interval is the flush frequency in milliseconds.
 * printdesc is either graphlab::event_log::NUMBER, or graphlab::event_log::DESCRIPTION or
 * graphlab::event_log::RATE_BAR
 * 
 * ADD_EVENT_TYPE(name, id, desc)
 * Creates an event type with an integer ID, and a description.
 * Event types should be mostly consecutive since the internal
 * storage format is a array. Valid ID range is [0, 255]
 *
 * ACCUMULATE_EVENT(name, id, count)
 * Adds 'count' events of type "id"
 *
 * FLUSH_EVENT_LOG(name)
 * Forces a flush of the accumulated events to the provided output stream
 */
#ifdef USE_EVENT_LOG
#define DECLARE_EVENT_LOG(name) graphlab::event_log name;
#define INITIALIZE_EVENT_LOG(name, ostrm, flush_interval, printdesc)    \
                    name.initialize(ostrm, flush_interval, printdesc);
#define ADD_EVENT_TYPE(name, id, desc) name.add_event_type(id, desc);
#define ACCUMULATE_EVENT(name, id, count) name.accumulate_event(id, count);
#define FLUSH_EVENT_LOG(name) name.flush();
#define CLOSE_EVENT_LOG(name) name.close();

#else
#define DECLARE_EVENT_LOG(name) 
#define INITIALIZE_EVENT_LOG(name, ostrm, flush_interval, printdesc)
#define ADD_EVENT_TYPE(name, id, desc) 
#define ACCUMULATE_EVENT(name, id, count) 
#define FLUSH_EVENT_LOG(name) 
#define CLOSE_EVENT_LOG(name)
#endif

#endif