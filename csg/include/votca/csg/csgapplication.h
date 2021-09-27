/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_CSG_CSGAPPLICATION_H
#define VOTCA_CSG_CSGAPPLICATION_H

// VOTCA includes
#include <memory>
#include <votca/tools/application.h>
#include <votca/tools/mutex.h>
#include <votca/tools/thread.h>

// Local VOTCA includes
#include "cgobserver.h"
#include "topology.h"
#include "topologymap.h"
#include "trajectoryreader.h"

namespace votca {
namespace csg {

class CsgApplication : public tools::Application {
 public:
  CsgApplication() = default;
  ~CsgApplication() override = default;

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

  virtual bool EvaluateTopology(Topology *, Topology * = nullptr) {
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

    /// \brief overload with the actual computation
    virtual void EvalConfiguration(Topology *top,
                                   Topology *top_ref = nullptr) = 0;

    /// \brief returns worker id
    Index getId() { return id_; }

   protected:
    CsgApplication *app_ = nullptr;
    Topology top_, top_cg_;
    std::unique_ptr<TopologyMap> map_;
    Index id_ = -1;

    void Run(void) override;

    void setApplication(CsgApplication *app) { app_ = app; }

    void setId(Index id) { id_ = id; }

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
  virtual std::unique_ptr<Worker> ForkWorker(void);

  /**
   * User is required to overload MergeWorker and merge data from each worker.
   * @param worker
   */
  virtual void MergeWorker(Worker *worker);

 protected:
  std::list<CGObserver *> observers_;
  bool do_mapping_;
  std::vector<std::unique_ptr<Worker>> myWorkers_;
  Index nframes_;
  bool is_first_frame_;
  Index nthreads_;
  tools::Mutex nframesMutex_;
  tools::Mutex traj_readerMutex_;

  /// \brief stores Mutexes used to impose order for input
  std::vector<std::unique_ptr<tools::Mutex>> threadsMutexesIn_;
  /// \brief stores Mutexes used to impose order for output
  std::vector<std::unique_ptr<tools::Mutex>> threadsMutexesOut_;
  std::unique_ptr<TrajectoryReader> traj_reader_;
};

inline void CsgApplication::AddObserver(CGObserver *observer) {
  observers_.push_back(observer);
}

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_CSGAPPLICATION_H
