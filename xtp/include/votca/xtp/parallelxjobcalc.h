/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// VOTCA includes
#include <votca/tools/mutex.h>

// Local VOTCA includes
#include "job.h"
#include "jobcalculator.h"
#include "progressobserver.h"
#include "qmthread.h"

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
  using Job = typename JobContainer::value_type;
  using Result = typename Job::JobResult;

  ParallelXJobCalc() = default;
  ~ParallelXJobCalc() override = default;

  std::string Identify() const override = 0;

  void ParseOptions(const tools::Property &opt) final {
    ParseCommonOptions(opt);
    ParseSpecificOptions(opt);
  }

  bool Evaluate(const Topology &top) final;
  virtual void CustomizeLogger(QMThread &thread);
  virtual Result EvalJob(const Topology &top, Job &job, QMThread &thread) = 0;

  void LockCout() { coutMutex_.Lock(); }
  void UnlockCout() { coutMutex_.Unlock(); }
  void LockLog() { logMutex_.Lock(); }
  void UnlockLog() { logMutex_.Unlock(); }

  // ======================================== //
  // XJOB OPERATOR (THREAD)                   //
  // ======================================== //

  class JobOperator : public QMThread {
   public:
    JobOperator(Index id, const Topology &top,
                ParallelXJobCalc<JobContainer> &master, Index openmp_threads)
        : top_(top), master_(master), openmp_threads_(openmp_threads) {
      setId(id);
    }  // comes from baseclass so Id cannot be in initializer list
    ~JobOperator() override = default;

    void Run() override;

   private:
    const Topology &top_;
    ParallelXJobCalc<JobContainer> &master_;
    Index openmp_threads_ = 1;
  };

 protected:
  virtual void ParseSpecificOptions(const tools::Property &options) = 0;

  // set the basis sets and functional in DFT package
  tools::Property UpdateDFTOptions(const tools::Property &options);
  // set the basis sets and functional in GWBSE
  tools::Property UpdateGWBSEOptions(const tools::Property &options);

  JobContainer XJobs_;
  tools::Mutex coutMutex_;
  tools::Mutex logMutex_;
  std::string mapfile_ = "";
  std::string jobfile_ = "";

 private:
  void ParseCommonOptions(const tools::Property &options);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_PARALLELXJOBCALC_H
