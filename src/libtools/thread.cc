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

#include "../../include/votca/tools/thread.h"
#include "../../include/votca/tools/lexical_cast.h"
#include <stdexcept>

using namespace std;

namespace votca {
namespace tools {

void *runwrapper(void *arg) {
  Thread *thread = (Thread *)(arg);
  thread->Run();
  pthread_exit(nullptr);
  return nullptr;
}

Thread::Thread() {}

Thread::~Thread() {}

void Thread::Start() {

  pthread_attr_init(&_attr);
  /*
   * according to the POSIX standard, threads are created in the joinable state
   * by default
   * however, some platforms do not obey
   *
   * explicitly create the thread in the joinable state
   * the main program can then wait for the thread by calling
   * pthread_join(thread)
   *
   */
  pthread_attr_setdetachstate(&_attr, PTHREAD_CREATE_JOINABLE);
  _finished = false;

  int rc =
      pthread_create(&_thread, &_attr, runwrapper, static_cast<void *>(this));
  if (rc) {
    throw std::runtime_error("ERROR; return code from pthread_create() is " +
                             boost::lexical_cast<std::string>(rc));
  }
}

void Thread::WaitDone() {
  void *status;
  int rc = pthread_join(_thread, &status);
  if (rc) {
    throw std::runtime_error("ERROR; return code from pthread_join() is " +
                             boost::lexical_cast<std::string>(rc));
  }
  _finished = true;
  pthread_attr_destroy(&_attr);
}

bool Thread::IsFinished() const { return _finished; }
}  // namespace tools
}  // namespace votca
