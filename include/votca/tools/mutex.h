/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef MUTEX_H
#define	MUTEX_H
#include <pthread.h>

namespace votca {
    namespace tools {

/**
         \brief Convenient class for Mutexes

         * Class allows to create, lock and unlock mutexes. Destroying is handled
         * by the destructor.

         */
        class Mutex {
        public:
            Mutex();
            ~Mutex();


            void Lock();
            void Unlock();

        private:
            pthread_mutex_t _mutexVar;
        };

    }
}

#endif	/* MUTEX_H */

