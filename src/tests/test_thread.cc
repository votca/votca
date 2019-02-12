/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE thread_test
#include "../../include/votca/tools/thread.h"
#include <boost/test/unit_test.hpp>
#include <cassert>
#include <exception>
#include <memory>
#include <vector>

using namespace std;
using namespace votca::tools;

class ThreadTest : public Thread {
 private:
  long factorial_;
  long sum_;

 public:
  long getFactorial() const { return sum_; }
  void setFactorial(long factorial) { factorial_ = factorial; }

 protected:
  void Run(void) {
    sum_ = 1;
    for (long val = 2; val <= factorial_; ++val) {
      sum_ *= val;
    }
  }
};

BOOST_AUTO_TEST_SUITE(thread_test)

BOOST_AUTO_TEST_CASE(constructors_test) { ThreadTest threadtest; }

BOOST_AUTO_TEST_CASE(start_to_finish_test) {
  ThreadTest threadtest;
  threadtest.setFactorial(3);
  threadtest.Start();
  threadtest.WaitDone();
  assert(threadtest.IsFinished());
  assert(threadtest.getFactorial() == 6);
}

BOOST_AUTO_TEST_CASE(multiple_start_to_finish_test) {

  vector<shared_ptr<ThreadTest>> threads;

  int numberThreads = 6;

  // Placing threads in a vector and initializing them
  for (int count = 0; count < numberThreads; ++count) {
    int factorial_value = count + 1;
    auto threadtest = make_shared<ThreadTest>();
    threadtest->setFactorial(factorial_value);
    threads.push_back(threadtest);
  }

  // Running the threads
  for (int count = 0; count < numberThreads; ++count) {
    threads.at(count)->Start();
  }

  // Wait until each of the threads is done
  for (int count = (numberThreads - 1); count >= 0; --count) {
    threads.at(count)->WaitDone();
  }

  // Check the factorials
  assert(threads.at(0)->getFactorial() == 1);
  assert(threads.at(1)->getFactorial() == 2);
  assert(threads.at(2)->getFactorial() == 6);
  assert(threads.at(3)->getFactorial() == 24);
  assert(threads.at(4)->getFactorial() == 120);
  assert(threads.at(5)->getFactorial() == 720);
}

BOOST_AUTO_TEST_SUITE_END()
