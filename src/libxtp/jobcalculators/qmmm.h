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

#pragma once
#ifndef VOTCA_XTP_QMMM_H
#define VOTCA_XTP_QMMM_H

// Local VOTCA includes
#include "votca/xtp/parallelxjobcalc.h"

namespace votca {
namespace xtp {

/**
 * \brief QM/MM with different regions around
 *
 * Calculates properties of different regions inside a multiregion calculation
 *
 * Callname: qmmm
 */

class QMMM final : public ParallelXJobCalc<std::vector<Job> > {
 public:
  std::string Identify() const { return "qmmm"; }
  Job::JobResult EvalJob(const Topology& top, Job& job, QMThread& Thread);
  void WriteJobFile(const Topology& top);
  void ReadJobFile(Topology& top);

 protected:
  void ParseSpecificOptions(const tools::Property& user_options);

 private:
  bool hasQMRegion() const;
  Job createJob(const Segment& seg, const QMState& state, Index jobid) const;
  std::string getFirstRegionName() const;

  std::pair<std::string, tools::Property> regions_def_;

  Index max_iterations_;
  bool print_regions_pdb_ = false;
  bool use_gs_for_ex_ = false;
  std::vector<QMState> states_;
  std::string which_segments_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_QMMM_H
