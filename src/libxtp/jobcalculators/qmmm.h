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

#pragma once
#ifndef VOTCA_XTP_QMMM_H
#define VOTCA_XTP_QMMM_H

#include <votca/xtp/parallelxjobcalc.h>

namespace votca {
namespace xtp {

/**
 * \brief QM/MM with different regions around
 *
 * Calculates properties of different regions inside a multiregion calculation
 *
 * Callname: qmmm
 */

class QMMM : public ParallelXJobCalc<std::vector<Job> > {
 public:
  void Initialize(tools::Property& options);
  std::string Identify() { return "qmmm"; }
  Job::JobResult EvalJob(Topology& top, Job& job, QMThread& Thread);
  void WriteJobFile(Topology& top);
  void ReadJobFile(Topology& top);

 private:
  tools::Property _regions_def;

  int _max_iterations = 100;
  bool _print_regions_pdb = false;
};

}  // namespace xtp
}  // namespace votca
#endif
