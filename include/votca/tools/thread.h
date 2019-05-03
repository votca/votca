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

#pragma once
#ifndef VOTCA_TOOLS_THREAD_H
#define VOTCA_TOOLS_THREAD_H

#include <pthread.h>

namespace votca {
namespace tools {

/**
 * \brief Framework for threaded execution
 *
 * The Thread class handles the threaded execution of VOTCA code.
 * In its current state, it is based on POSIX threads.
 * It mainly contains a wrapper, a start and a wait function.
 * This class should not be touched by the user.
 **/
class Thread {
 public:
  Thread();
  virtual ~Thread();

  /**
   * \brief Starts and runs a thread
   *
   * This method is responsible for creating and running the thread.
   **/
  void Start();

  /**
   * \brief WaitDone() will not exit until thread ends computation.
   */
  void WaitDone();

  /**
   * \brief Checks to see if the thread is done
   *
   * \return True if thread is done, False otherwise
   **/
  bool IsFinished() const;

 protected:
  /**
   * \brief Run() executes the actual code.
   *
   * It is left blank here and must be overloaded, for the thread to do
   * any work. This method should never be called by the user, it is
   * called internally by the Start method.
   */
  virtual void Run(void) = 0;
  friend void *runwrapper(void *arg);

 private:
  pthread_attr_t _attr;
  pthread_t _thread;
  bool _finished;
};
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_THREAD_H
