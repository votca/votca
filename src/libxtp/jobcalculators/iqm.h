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
#ifndef VOTCA_XTP_IQM_H
#define VOTCA_XTP_IQM_H

#include <votca/tools/property.h>

#include <boost/filesystem.hpp>
#include <sys/stat.h>
#include <votca/xtp/bsecoupling.h>
#include <votca/xtp/dftcoupling.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/parallelxjobcalc.h>

namespace votca {
namespace xtp {

/**
 * \brief DFT & GWBSE-based coupling elements
 *
 * Evaluates DFT & GWBSE-based coupling elements for all conjugated
 * segments from the neighbor list. Requires molecular orbitals of two monomers
 * and a dimer in GAUSSIAN, NWChem, or ORCAformat.
 *
 * Callname: iqm
 */

class IQM : public ParallelXJobCalc<std::vector<Job> > {
 public:
  void Initialize(tools::Property& options) override;
  std::string Identify() override { return "iqm"; }
  Job::JobResult EvalJob(const Topology& top, Job& job,
                         QMThread& Thread) override;
  void WriteJobFile(const Topology& top) override;
  void ReadJobFile(Topology& top) override;

 private:
  double GetBSECouplingFromProp(tools::Property& bseprop, const QMState& stateA,
                                const QMState& stateB);
  double GetDFTCouplingFromProp(tools::Property& dftprop, Index stateA,
                                Index stateB);
  void SetJobToFailed(Job::JobResult& jres, Logger& pLog,
                      const std::string& errormessage);
  void WriteLoggerToFile(const std::string& logfile, Logger& logger);
  void addLinkers(std::vector<const Segment*>& segments, const Topology& top);
  bool isLinker(const std::string& name);
  void ParseOptionsXML(tools::Property& opt);
  std::map<std::string, QMState> FillParseMaps(const std::string& Mapstring);

  QMState GetElementFromMap(const std::map<std::string, QMState>& elementmap,
                            const std::string& elementname) const;

  tools::Property _dftpackage_options;
  tools::Property _gwbse_options;
  tools::Property _bsecoupling_options;
  tools::Property _dftcoupling_options;

  // what to do
  bool _do_dft_input = false;
  bool _do_dft_run = false;
  bool _do_dft_parse = false;
  bool _do_dftcoupling = false;
  bool _do_gwbse = false;
  bool _do_bsecoupling = false;

  std::map<std::string, QMState> _linkers;

  // what to write in the storage
  bool _store_dft = false;
  bool _store_gw = false;

  // parsing options
  std::map<std::string, QMState> _singlet_levels;
  std::map<std::string, QMState> _triplet_levels;

  std::map<std::string, QMState> _hole_levels;
  std::map<std::string, QMState> _electron_levels;
};

}  // namespace xtp
}  // namespace votca
#endif
