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
#ifndef VOTCA_XTP_DFTGWBSE_H
#define VOTCA_XTP_DFTGWBSE_H

// Standard includes
#include <cstdio>

// Local VOTCA includes
#include "votca/xtp/dftengine.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {

class DftGwBse final : public QMTool {
 public:
  DftGwBse() = default;

  ~DftGwBse() = default;

  std::string Identify() const { return "dftgwbse"; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Run();

 private:
  std::string guess_file_;
  // Separate from guess_file_ above: guess_file_ loads BOTH geometry
  // and MOs from the same .orb checkpoint (via Orbitals::ReadFromCpt),
  // which cannot represent "new geometry, but warm-started MOs from a
  // previous, different geometry" -- exactly what a geometry
  // optimization needs between steps. moguess_file_ loads ONLY the
  // MOs from a separate .orb file, leaving the geometry itself to come
  // from xyzfile_ as usual. See this class's own Run() for how the
  // two are combined.
  std::string moguess_file_;
  std::string mpsfile_;
  // QMState (e.g. "s1") whose analytic forces are computed after the BSE
  std::string excited_state_forces_;
  bool excited_state_forces_qp_response_ = true;
  // orb file of an earlier geometry: follow its excited_state_forces state
  std::string excited_state_forces_follow_;
  void ComputeExcitedStateForces(Orbitals& orbitals);
  // QP levels whose energy gradients are computed, e.g. homo-2:lumo+3
  std::string qp_gradients_;
  void ComputeQPGradients(Orbitals& orbitals);
  // tasks = gradients: gradients from the orb file of an earlier run
  bool gradients_only_ = false;
  void CheckStoredRun(const Orbitals& orbitals) const;
  void SetupQPResponse(DFTEngine::BSEGradientInput& input,
                       const Orbitals& orbitals);
  static std::vector<Index> ParseLevels(const std::string& text, Index homo);
  QMState FollowState(const Orbitals& orbitals);

  std::string xyzfile_;
  std::string xml_output_;    // .xml output
  std::string archive_file_;  // .orb file to parse to

  tools::Property package_options_;
  tools::Property gwbseengine_options_;
  tools::Property geoopt_options_;

  Logger log_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTGWBSE_H
