/*
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_CSGAPPLICATION_H
#define	__VOTCA_CSGAPPLICATION_H

#include <votca/tools/application.h>
#include "topology.h"
#include "topologymap.h"
#include "cgobserver.h"
#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>
#include "trajectoryreader.h"

namespace votca {
    namespace csg {
        using namespace votca::tools;

        class CsgApplication
        : public Application {
        public:
            CsgApplication();
            ~CsgApplication();

            void Initialize();
            bool EvaluateOptions();

            void Run(void);
            void RunThreaded(void);

            void ShowHelpText(std::ostream &out);

            /// \brief overload and return true to enable mapping command line options

            virtual bool DoMapping(void) {
                return false;
            }
            /// \brief if DoMapping is true, will by default require mapping or not

            virtual bool DoMappingDefault(void) {
                return true;
            }
            /// \brief overload and return true to enable trajectory command line options

            virtual bool DoTrajectory(void) {
                return false;
            }

            /* \brief overload and return true to enable threaded calculations */
            virtual bool DoThreaded(void) {
                return false;
            }

            /* \brief overload and return false to disable synchronized (while threaded) calculations */
            virtual bool SynchronizeThreads(void) {
                if (DoThreaded())
                    return true;
                else
                    return false;
            }

            /// \brief called after topology was loaded

            virtual bool EvaluateTopology(Topology *top, Topology *top_ref = 0) {
                return true;
            }

            void AddObserver(CGObserver *observer);

            /// \brief called before the first frame
            virtual void BeginEvaluate(Topology *top, Topology *top_ref = 0);
            /// \brief called after the last frame
            virtual void EndEvaluate();
            // \brief called for each frame which is mapped
            virtual void EvalConfiguration(Topology *top, Topology *top_ref = 0);


            // thread related stuff follows

            class Worker : public Thread {
            public:

                Worker();
                ~Worker();

                virtual void EvalConfiguration(Topology *top, Topology *top_ref = 0) = 0;

                int getId() {
                    return _id;
                }

            protected:
                CsgApplication *_app;
                Topology _top, _top_cg;
                TopologyMap * _map;
                int _id;

                void Run(void);

                void setApplication(CsgApplication *app) {
                    _app = app;
                }

                void setId(int id) {
                    _id = id;
                }

                friend class CsgApplication;
            };

            /**
             *  TODO comment!
             * @param worker
             * @return 
             */
            bool ProcessData(Worker * worker);

            /**
             * TODO comment
             * @return
             */
            virtual Worker *ForkWorker(void);

            /**
             * TODO comment
             * @param worker
             */
            virtual void MergeWorker(Worker *worker);

        protected:
            list<CGObserver *> _observers;
            bool _do_mapping;
            std::vector<Worker*> _myWorkers;
            int _nframes;
            bool _is_first_frame;
            int _nthreads;
            Mutex _nframesMutex;
            Mutex _traj_readerMutex;
            std::vector<Mutex*> _threadsMutexesIn;
            std::vector<Mutex*> _threadsMutexesOut;
            TrajectoryReader * _traj_reader;
        };

        inline void CsgApplication::AddObserver(CGObserver *observer) {
            _observers.push_back(observer);
        }

        inline CsgApplication::Worker::Worker()
        : _map(NULL), _id(-1), _app(NULL) {
        }

    }
}

#endif	/* __VOTCA_CSGAPPLICATION_H */

