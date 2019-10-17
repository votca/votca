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

#ifndef _VOTCA_CSG_APPLICATION_H
#define _VOTCA_CSG_APPLICATION_H

#include "cgobserver.h"
#include "topology.h"
#include "topologymap.h"
#include "trajectoryreader.h"
#include <votca/tools/application.h>
#include <votca/tools/mutex.h>
#include <votca/tools/thread.h>

namespace votca {
namespace csg {

class CsgApplication : public tools::Application {
 public:
  CsgApplication() = default;
  ;
  ~CsgApplication() override = default;
  ;

  void Initialize() override;
  bool EvaluateOptions() override;

  void Run(void) override;

  void ShowHelpText(std::ostream &out) override;

  /// \brief overload and return true to enable mapping command line options

  virtual bool DoMapping(void) { return false; }
  /// \brief if DoMapping is true, will by default require mapping or not

  virtual bool DoMappingDefault(void) { return true; }
  /// \brief overload and return true to enable trajectory command line options

  virtual bool DoTrajectory(void) { return false; }

  /* \brief overload and return true to enable threaded calculations */
  virtual bool DoThreaded(void) { return false; }

  /* \brief overload and return false to disable synchronized (while threaded)
   * calculations */
  virtual bool SynchronizeThreads(void) {
    if (DoThreaded()) {
      return true;
    } else {
      return false;
    }
  }

  /// \brief if topology is always needed
  virtual bool NeedsTopology(void) { return true; }

  /// \brief called after topology was loaded

  virtual bool EvaluateTopology(Topology *top, Topology *top_ref = nullptr) {
    return true;
  }

  void AddObserver(CGObserver *observer);

  /// \brief called before the first frame
  virtual void BeginEvaluate(Topology *top, Topology *top_ref = nullptr);
  /// \brief called after the last frame
  virtual void EndEvaluate();
  // \brief called for each frame which is mapped
  virtual void EvalConfiguration(Topology *top, Topology *top_ref = nullptr);

  // thread related stuff follows

  /**
   \brief Worker, derived from Thread, does the work.
   *
   * Worker holds the information about the current frame, either in its
   * own copy (e.g. Topology), or, by reference, from the parent CsgApplication.
   * The computation is shifted from Run() into EvalConfiguration. The
   * user is required to overload ForkWorker and Mergeworker and thereby
   * define the initialization and merging of workers. By default, workers
   * will be executed in correct order according to the frames. Also,
   * output will follow the same order.
   * Mutexes handle the locking of input/output and are also used to impose
   * the correct order of frames for in/output.
   *
   */
  class Worker : public tools::Thread {
   public:
    Worker() = default;
    ;
    ~Worker() override;

    /// \brief overload with the actual computation
    virtual void EvalConfiguration(Topology *top,
                                   Topology *top_ref = nullptr) = 0;

    /// \brief returns worker id
    int getId() { return _id; }

   protected:
    CsgApplication *_app = nullptr;
    Topology _top, _top_cg;
    TopologyMap *_map = nullptr;
    int _id = -1;

    void Run(void) override;

    void setApplication(CsgApplication *app) { _app = app; }

    void setId(int id) { _id = id; }

    friend class CsgApplication;
  };

  /**
   * \brief Gets frames from TrajectoryReader in an ordered way and, if
   * successful, calls Worker::EvalConfiguration for that frame.
   *
   * @param worker
   * @return True if frames left for calculation, else False
   */
  bool ProcessData(Worker *worker);

  /**
   *
   * User is required to overload ForkWorker and initialize workers.
   * @return worker
   */
  virtual Worker *ForkWorker(void);

  /**
   * User is required to overload MergeWorker and merge data from each worker.
   * @param worker
   */
  virtual void MergeWorker(Worker *worker);

 protected:
  std::list<CGObserver *> _observers;
  bool _do_mapping;
  std::vector<Worker *> _myWorkers;
  int _nframes;
  bool _is_first_frame;
  int _nthreads;
  tools::Mutex _nframesMutex;
  tools::Mutex _traj_readerMutex;

  /// \brief stores Mutexes used to impose order for input
  std::vector<tools::Mutex *> _threadsMutexesIn;
  /// \brief stores Mutexes used to impose order for output
  std::vector<tools::Mutex *> _threadsMutexesOut;
  TrajectoryReader *_traj_reader;
};

inline void CsgApplication::AddObserver(CGObserver *observer) {
  _observers.push_back(observer);
}

}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_APPLICATION_H */
