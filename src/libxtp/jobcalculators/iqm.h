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

#ifndef _CALC_XTP_IQM_H
#define _CALC_XTP_IQM_H

#include <votca/tools/property.h>

#include <boost/filesystem.hpp>
#include <sys/stat.h>
#include <votca/ctp/parallelxjobcalc.h>
#include <votca/xtp/bsecoupling.h>
#include <votca/xtp/dftcoupling.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/orbitals.h>

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

class IQM : public ctp::ParallelXJobCalc<vector<ctp::Job*>, ctp::Job*,
                                         ctp::Job::JobResult> {
 public:
  void Initialize(tools::Property* options);
  std::string Identify() { return "iqm"; }
  ctp::Job::JobResult EvalJob(ctp::Topology* top, ctp::Job* job,
                              ctp::QMThread* Thread);
  void WriteJobFile(ctp::Topology* top);
  void ReadJobFile(ctp::Topology* top);

 private:
  double GetBSECouplingFromProp(tools::Property& bseprop, const QMState& stateA,
                                const QMState& stateB);
  double GetDFTCouplingFromProp(tools::Property& dftprop, int stateA,
                                int stateB);
  void SetJobToFailed(ctp::Job::JobResult& jres, ctp::Logger* pLog,
                      const std::string& errormessage);
  void WriteLoggerToFile(const std::string& logfile, ctp::Logger& logger);
  void addLinkers(std::vector<ctp::Segment*>& segments, ctp::Topology* top);
  bool isLinker(const std::string& name);
  void WriteCoordinatesToOrbitalsPBC(ctp::QMPair& pair, Orbitals& orbitals);
  void ParseOptionsXML(tools::Property& opt);
  std::map<std::string, QMState> FillParseMaps(const string& Mapstring);

  QMState GetElementFromMap(const std::map<std::string, QMState>& elementmap,
                            const std::string& elementname) const;

  std::string _package;
  Property _dftpackage_options;
  Property _gwbse_options;
  Property _bsecoupling_options;
  Property _dftcoupling_options;

  // what to do
  bool _do_dft_input = false;
  bool _do_dft_run = false;
  bool _do_dft_parse = false;
  bool _do_dftcoupling = false;
  bool _do_gwbse = false;
  bool _do_bsecoupling = false;

  std::vector<std::string> _linker_names;

  // what to write in the storage
  bool _store_dft = false;
  bool _store_singlets = false;
  bool _store_triplets = false;
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
