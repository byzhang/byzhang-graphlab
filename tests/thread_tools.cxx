#include <iostream>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/parallel/thread_pool.hpp>
#include <graphlab/parallel/thread_flip_flop.hpp>
#include <graphlab/logger/assertions.hpp>
#include <boost/bind.hpp>

using namespace graphlab;

atomic<int> testval;

void test_inc() {
  usleep(100000);
  testval.inc();
}

void test_dec() {
  usleep(100000);
  testval.dec();
}



void thread_assert_false() {
  ASSERT_TRUE(false);
}


void test_group_exception_forwarding(){
  std::cout << "\n";
  std::cout << "----------------------------------------------------------------\n";
  std::cout << "This test will print a  large number of assertional failures\n";
  std::cout << "and back traces. This is intentional as we are testing the\n" ;
  std::cout << "exception forwarding scheme\n";
  std::cout << "----------------------------------------------------------------\n";
  std::cout << std::endl;

  thread_group group;

  
  thread thr3;
  thr3.launch(thread_assert_false);
  try {
    thr3.join();
  }
  catch(const char* c) {
    std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
  }
  
  
  for (size_t i = 0;i < 10; ++i) {
    group.launch(thread_assert_false);
  }
  
  size_t numcaught = 0;
  while (group.running_threads() > 0) {
    try {
      group.join();
    }
    catch (const char* c){
      std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
      numcaught++;
    }
  }
  std::cout << "Caught " << numcaught << " exceptions!" << std::endl;
  TS_ASSERT_EQUALS(numcaught, (size_t)10);
}

void test_pool(){
  testval.value = 0;
  thread_pool pool(4);
  for (size_t j = 0;j < 10; ++j) {
    for (size_t i = 0;i < 10; ++i) {
      pool.launch(test_inc);
    }
    for (size_t i = 0;i < 10; ++i) {
      pool.launch(test_dec);
    }
    pool.set_cpu_affinity(j % 2);
  }
  
  pool.join();
  TS_ASSERT_EQUALS(testval.value, 0);
}

void test_pool_exception_forwarding(){
  std::cout << "\n";
  std::cout << "----------------------------------------------------------------\n";
  std::cout << "This test will print a  large number of assertional failures\n";
  std::cout << "and back traces. This is intentional as we are testing the\n" ;
  std::cout << "exception forwarding scheme\n";
  std::cout << "----------------------------------------------------------------\n";
  std::cout << std::endl;
  thread_pool pool(10);

  
  thread thr3;
  thr3.launch(thread_assert_false);
  try {
    thr3.join();
  }
  catch(const char* c) {
    std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
  }
  
  
  for (size_t i = 0;i < 10; ++i) {
    pool.launch(thread_assert_false);
    if (i == 50) {
      pool.set_cpu_affinity(true);
    }
  }
  
  size_t numcaught = 0;
  while (1) {
    try {
      pool.join();
      break;
    }
    catch (const char* c){
      std::cout << "Exception " << c << " forwarded successfully!" << std::endl;
      numcaught++;
    }
  }
  std::cout << "Caught " << numcaught << " exceptions!" << std::endl;
  TS_ASSERT_EQUALS(numcaught, (size_t)10);
}



const int num_flip_flop_threads = 4;
barrier f0_barrier(num_flip_flop_threads);
barrier f1_barrier(num_flip_flop_threads);
atomic<size_t> f0;
atomic<size_t> f1;
thread_flip_flop flipflop(num_flip_flop_threads,
                          num_flip_flop_threads);

void flip_flop_0() {
  for (size_t i = 0; i < 1000; ++i) {
    for (size_t j = 0;j < 10; ++j) {
      f0.inc();
    }
    f0_barrier.wait();
    ASSERT_EQ(f0.value - num_flip_flop_threads * 10, f1.value);
    flipflop.wait(0);
  }
  flipflop.stop_blocking();
}

void flip_flop_1() {
  for (size_t i = 0; i < 1000; ++i) {
    flipflop.wait(1);
    for (size_t j = 0;j < 10; ++j) {
      f1.inc();
    }
  }
  flipflop.stop_blocking();
}

class ThreadToolsTestSuite : public CxxTest::TestSuite {
public:
   void test_thread_group_exception(void) {
    test_group_exception_forwarding();
   }

   void test_thread_pool(void) {
    test_pool();
   }
   
   void test_thread_pool_exception(void) {
     test_pool_exception_forwarding();
   }

};
