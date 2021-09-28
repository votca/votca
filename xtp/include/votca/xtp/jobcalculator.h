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
#ifndef VOTCA_XTP_JOBCALCULATOR_H
#define VOTCA_XTP_JOBCALCULATOR_H

// VOTCA includes
#include <votca/tools/calculator.h>

// Local VOTCA includes
#include "job.h"
#include "progressobserver.h"
#include "topology.h"

namespace votca {
namespace xtp {

class Topology;

class JobCalculator : public tools::Calculator {
 public:
  JobCalculator() = default;
  ~JobCalculator() override = default;

  std::string Identify() const override = 0;
  std::string Package() const final { return "xtp"; }
  bool EvaluateFrame(const Topology &top) { return Evaluate(top); }

  void Initialize(const tools::Property &opt) final { ParseOptions(opt); }
  virtual void WriteJobFile(const Topology &top) = 0;
  virtual void ReadJobFile(Topology &top) = 0;

  void setOpenMPThreads(Index ompthreads) { openmp_threads_ = ompthreads; }
  void setProgObserver(ProgObserver<std::vector<Job> > *obs) { progObs_ = obs; }

 protected:
  virtual void ParseOptions(const tools::Property &opt) = 0;
  virtual bool Evaluate(const Topology &top) = 0;

  Index openmp_threads_;
  ProgObserver<std::vector<Job> > *progObs_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_JOBCALCULATOR_H
