/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#ifndef VOTCA_XTP_PODCOUPLINGTOOL_H
#define VOTCA_XTP_PODCOUPLINGTOOL_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/podcoupling.h"
#include "votca/xtp/qmtool.h"
#include <votca/tools/types.h>

namespace votca {
namespace xtp {

// Standalone tool wrapping PODCoupling -- direct, low-risk first piece
// of the calculator integration discussed with the user, matching the
// established QMTool pattern (Coupling, Diabatization) rather than
// either of those tools' own, structurally different input shapes:
// unlike Coupling (three separate quantum-chemistry outputs: two
// monomers + a dimer) or Diabatization (one or two ADIABATIC states,
// with a state_idx each), this tool consumes a SINGLE, already-
// converged .orb file for the intact supermolecule, plus two explicit
// atom-index fragment definitions -- matching PODCoupling's own,
// genuinely different input shape directly rather than forcing it
// into either existing tool's own convention.
class PodCouplingTool final : public QMTool {
 public:
  PodCouplingTool() = default;
  ~PodCouplingTool() = default;

  std::string Identify() const { return "podcoupling"; }

 protected:
  void ParseOptions(const tools::Property& user_options) final;
  bool Run() final;

 private:
  std::string orb_file_;
  Logger log_;

  std::vector<Index> fragment_A_atoms_;
  std::vector<Index> fragment_B_atoms_;

  // Range of orbitals to cover for each fragment, matching
  // DFTcoupling's own numberofstatesA_/numberofstatesB_ convention
  // exactly (XML option names levA/levB, same as dftcoupling.cc's own
  // options.get("levA")/("levB")) -- per direct agreement with the
  // user to mirror DFTcoupling's own established behavior here: BOTH
  // occupied (hole-transport) and virtual (electron-transport)
  // orbitals are covered together, in a single call, rather than this
  // tool's own, earlier, single-orbital-pair-only interface. N=1
  // covers just each fragment's own {HOMO, LUMO}.
  Index numberofstatesA_;
  Index numberofstatesB_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_PODCOUPLINGTOOL_H
