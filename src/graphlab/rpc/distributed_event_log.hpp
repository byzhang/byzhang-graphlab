#ifndef DISTRIBUTED_EVENT_LOG_HPP
#define DISTRIBUTED_EVENT_LOG_HPP
#include <iostream>
#include <boost/bind.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/atomic.hpp>
#include <graphlab/util/timer.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
namespace graphlab {

#define EVENT_MAX_COUNTERS 256

class dist_event_log{
 public:
  enum event_print_type{
    NUMBER,
    DESCRIPTION,
    RATE_BAR
  };
 private:
  std::ostream* out;  // output target
  volatile size_t flush_interval; // flush frequency (ms)
  event_print_type print_method;  // print method
  volatile bool finished; // set on destruction
  
  mutex m;          
  conditional cond;
  
  thread printing_thread;

  // Local Counters
  std::string descriptions[EVENT_MAX_COUNTERS];
  atomic<size_t> counters[EVENT_MAX_COUNTERS];

  // Global counters on proc 0
  std::vector<atomic<size_t> > globalcounters[EVENT_MAX_COUNTERS];
  size_t maxcounter[EVENT_MAX_COUNTERS];
  size_t maxproc_counter[EVENT_MAX_COUNTERS];  // maximum of per proc maximums

  timer ti;
  double prevtime;  // last time flush() was called
  bool hasevents;   // whether there are logging events
  size_t noeventctr; // how many times flush() was called with no events
                     // zeroed when the events are next observed.

  size_t max_desc_length; // maximum descriptor length
  fixed_dense_bitset<EVENT_MAX_COUNTERS> hascounter;

  dc_dist_object<dist_event_log> *rmi;
 public:
  inline dist_event_log():out(NULL),
                        flush_interval(0),
                        print_method(DESCRIPTION),
                        finished(false),
                        prevtime(0),
                        hasevents(false),
                        noeventctr(0),
                        max_desc_length(0),
                        rmi(NULL) {
    hascounter.clear();
    for (size_t i = 0;i < EVENT_MAX_COUNTERS; ++i) {
      maxcounter[i] = 0;
      maxproc_counter[i] = 0;
    }
    printing_thread.launch(boost::bind(&dist_event_log::thread_loop, this));
  }
  
  void initialize(distributed_control& dc,
                  std::ostream &ostrm,
                  size_t flush_interval_ms,
                  event_print_type event_print);

  void close();
  
  void thread_loop();

  void add_event_type(unsigned char eventid, std::string description);

  void accumulate_event_aggregator(procid_t proc,
                                   unsigned char eventid,
                                   size_t count);
  
  inline void accumulate_event(unsigned char eventid,
                               size_t count)  __attribute__((always_inline)) {
    counters[eventid].inc(count);
  }

  void print_log();
  void flush();
  ~dist_event_log();
};

}
/**
 * DECLARE_DIST_EVENT_LOG(name)
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
 * INITIALIZE_DIST_EVENT_LOG(name, dc, ostrm, flush_interval, printdesc)
 * ostrm is the output std::ostream object. A pointer of the stream
 * is taken so the stream should not be destroyed until the event log is closed
 * flush_interval is the flush frequency in milliseconds.
 * printdesc is either graphlab::event_log::NUMBER, or graphlab::event_log::DESCRIPTION or
 * graphlab::event_log::RATE_BAR. This must be called on all machines,
 * but only machine 0 will use the ostrm, and printdesc arguments
 * 
 * ADD_DIST_EVENT_TYPE(name, id, desc)
 * Creates an event type with an integer ID, and a description.
 * Event types should be mostly consecutive since the internal
 * storage format is a array. Valid ID range is [0, 255]
 *
 * ACCUMULATE_DIST_EVENT(name, id, count)
 * Adds 'count' events of type "id"
 *
 * FLUSH_DIST_EVENT_LOG(name)
 * Forces a flush of the accumulated events to the provided output stream
 */
#ifdef USE_EVENT_LOG
#define DECLARE_DIST_EVENT_LOG(name) graphlab::dist_event_log name;
#define INITIALIZE_DIST_EVENT_LOG(name, dc, ostrm, flush_interval, printdesc)    \
                    name.initialize(dc, ostrm, flush_interval, printdesc);
#define ADD_DIST_EVENT_TYPE(name, id, desc) name.add_event_type(id, desc);
#define ACCUMULATE_DIST_EVENT(name, id, count) name.accumulate_event(id, count);
#define FLUSH_DIST_EVENT_LOG(name) name.flush();
#define CLOSE_DIST_EVENT_LOG(name) name.close();

#else
#define DECLARE_DIST_EVENT_LOG(name) 
#define INITIALIZE_DIST_EVENT_LOG(name, dc, ostrm, flush_interval, printdesc)
#define ADD_DIST_EVENT_TYPE(name, id, desc) 
#define ACCUMULATE_DIST_EVENT(name, id, count) 
#define FLUSH_DIST_EVENT_LOG(name) 
#define CLOSE_DIST_EVENT_LOG(name)
#endif

#endif