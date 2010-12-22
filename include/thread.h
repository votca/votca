/*
 * Copyright 2010 The VOTCA Development Team (http://www.votca.org)
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

#ifndef THREAD_H
#define	THREAD_H

#include <pthread.h>

namespace votca {
    namespace tools {

        /**
            \brief Framework for threaded execution

         * The Thread class handles the threaded execution of VOTCA code.
         * In its current state, it is based on POSIX threads.
         * It mainly contains a wrapper, a start and a wait function.
         * This class should not be touched by the user.

         */
        class Thread {
        public:
            Thread();
            virtual ~Thread();

            void Start();

            /**
             * \brief Run() executes the actual code. Overload it.
             */
            virtual void Run(void) = 0;

            /**
             * \brief WaitDone() will not exit until thread ends computation.
             */
            void WaitDone();
            bool IsFinished() const;

        private:
            pthread_t _thread;
            bool _finished;

        };

    }
}

#endif	/* THREAD_H */

