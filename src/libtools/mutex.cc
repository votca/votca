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

#include <votca/tools/mutex.h>

namespace votca {
namespace tools {

Mutex::Mutex() { pthread_mutex_init(&_mutexVar, nullptr); }

Mutex::~Mutex() { pthread_mutex_destroy(&_mutexVar); }

void Mutex::Lock() { pthread_mutex_lock(&_mutexVar); }

void Mutex::Unlock() { pthread_mutex_unlock(&_mutexVar); }

}  // namespace tools
}  // namespace votca
