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
#ifndef VOTCA_XTP_IQM_H
#define VOTCA_XTP_IQM_H

// Third party includes
#include <boost/filesystem.hpp>
#include <sys/stat.h>

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/bsecoupling.h"
#include "votca/xtp/dftcoupling.h"
#include "votca/xtp/gwbse.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/parallelxjobcalc.h"

namespace votca {
namespace xtp {

/**
 * \brief DFT & GWBSE-based coupling elements
 *
 * Evaluates DFT & GWBSE-based coupling elements for all conjugated
 * segments from the neighbor list. Requires molecular orbitals of two monomers
 * and a dimer in ORCA format.
 *
 * Callname: iqm
 */

class IQM final : public ParallelXJobCalc<std::vector<Job> > {
 public:
  std::string Identify() { return "iqm"; }
  Job::JobResult EvalJob(const Topology& top, Job& job, QMThread& opThread);
  void WriteJobFile(const Topology& top);
  void ReadJobFile(Topology& top);

 protected:
  void ParseSpecificOptions(const tools::Property& user_options);

 private:
  double GetBSECouplingFromProp(const tools::Property& bseprop,
                                const QMState& stateA, const QMState& stateB);
  double GetDFTCouplingFromProp(const tools::Property& dftprop, Index stateA,
                                Index stateB);
  void SetJobToFailed(Job::JobResult& jres, Logger& pLog,
                      const std::string& errormessage);
  void WriteLoggerToFile(const std::string& logfile, Logger& logger);
  void addLinkers(std::vector<const Segment*>& segments, const Topology& top);
  bool isLinker(const std::string& name);
  std::map<std::string, QMState> FillParseMaps(const std::string& Mapstring);

  QMState GetElementFromMap(const std::map<std::string, QMState>& elementmap,
                            const std::string& elementname) const;

  tools::Property dftpackage_options_;
  tools::Property gwbse_options_;
  tools::Property bsecoupling_options_;
  tools::Property dftcoupling_options_;

  // what to do
  bool do_dft_input_ = false;
  bool do_dft_run_ = false;
  bool do_dft_parse_ = false;
  bool do_dftcoupling_ = false;
  bool do_gwbse_ = false;
  bool do_bsecoupling_ = false;

  std::map<std::string, QMState> linkers_;

  // what to write in the storage
  bool store_dft_ = false;
  bool store_gw_ = false;

  // parsing options
  std::map<std::string, QMState> singlet_levels_;
  std::map<std::string, QMState> triplet_levels_;

  std::map<std::string, QMState> hole_levels_;
  std::map<std::string, QMState> electron_levels_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_IQM_H
