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
#ifndef VOTCA_XTP_EQM_H
#define VOTCA_XTP_EQM_H

// Local VOTCA includes
#include "votca/xtp/gwbse.h"
#include "votca/xtp/parallelxjobcalc.h"
#include "votca/xtp/qmpackagefactory.h"
#include "votca/xtp/segment.h"

namespace votca {
namespace xtp {

/**
 * \brief Run DFT/GWBSE calculations
 *
 * Evaluates DFT and GWBSE for all molecules
 * Requires a first-principles package, i.e. ORCA
 *
 */

class EQM : public ParallelXJobCalc<std::vector<Job> > {
 public:
  std::string Identify() override { return "eqm"; }
  void Initialize(const tools::Property &user_options) override;
  Job::JobResult EvalJob(const Topology &top, Job &job,
                         QMThread &opThread) override;

  void CleanUp() { ; }
  void WriteJobFile(const Topology &top) override;
  void ReadJobFile(Topology &) override { return; }

 private:
  void WriteLoggerToFile(const std::string &logfile, Logger &logger);

  void SetJobToFailed(Job::JobResult &jres, Logger &pLog,
                      const std::string &errormessage);

  tools::Property _package_options;
  tools::Property _gwbse_options;
  tools::Property _esp_options;

  // what to do
  bool _do_dft_input = false;
  bool _do_dft_run = false;
  bool _do_dft_parse = false;
  bool _do_gwbse = false;
  bool _do_esp = false;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EQM_H
