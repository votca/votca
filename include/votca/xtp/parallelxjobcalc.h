/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
/// For an earlier history see ctp repo commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#pragma once
#ifndef VOTCA_XTP_PARALLELXJOBCALC_H
#define VOTCA_XTP_PARALLELXJOBCALC_H

#include <votca/tools/mutex.h>
#include <votca/xtp/job.h>
#include <votca/xtp/jobcalculator.h>
#include <votca/xtp/progressobserver.h>
#include <votca/xtp/qmthread.h>

/// PATHWAYS TO A NEW THREADED CALCULATOR
/// ... 1 Define 'JobContainer' (needs to define iterator), 'pJob' ( =
/// *iterator)
/// ... 2 Derive new calculator as ': public
/// ParallelXJobCalc<JobContainer,pJob>'
/// ... 3 Specialize XJOBS_FROM_TABLE< JobContainer, pJob> in xjob.cc
/// ... 4 Register new calculator (see end of parallelxjobcalc.cc)

/// REQUIRED MEMBERS FOR pJob
/// pJob::JobResult (struct)

namespace votca {
namespace xtp {

template <typename JobContainer>
class ParallelXJobCalc : public JobCalculator {

 public:
  class JobOperator;
  using Job = typename std::iterator_traits<
      typename JobContainer::iterator>::value_type;
  using Result = typename Job::JobResult;

  ParallelXJobCalc() = default;
  ~ParallelXJobCalc() override { ; };

  std::string Identify() override = 0;

  bool EvaluateFrame(const Topology &top) override;
  virtual void CustomizeLogger(QMThread &thread);
  virtual Result EvalJob(const Topology &top, Job &job, QMThread &thread) = 0;

  void LockCout() { _coutMutex.Lock(); }
  void UnlockCout() { _coutMutex.Unlock(); }
  void LockLog() { _logMutex.Lock(); }
  void UnlockLog() { _logMutex.Unlock(); }

  // ======================================== //
  // XJOB OPERATOR (THREAD)                   //
  // ======================================== //

  class JobOperator : public QMThread {
   public:
    JobOperator(Index id, const Topology &top,
                ParallelXJobCalc<JobContainer> &master, Index openmp_threads)
        : _top(top), _master(master), _openmp_threads(openmp_threads) {
      setId(id);
    }  // comes from baseclass so Id cannot be in initializer list
    ~JobOperator() override = default;

    void Run() override;

   private:
    const Topology &_top;
    ParallelXJobCalc<JobContainer> &_master;
    Index _openmp_threads = 1;
  };

 protected:
  void ParseCommonOptions(const tools::Property &options);

  JobContainer _XJobs;
  tools::Mutex _coutMutex;
  tools::Mutex _logMutex;
  std::string _mapfile = "";
  std::string _jobfile = "";
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_PARALLELXJOBCALC_H
