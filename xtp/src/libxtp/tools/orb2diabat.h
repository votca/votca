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
#ifndef VOTCA_XTP_ORB2DIABAT_H
#define VOTCA_XTP_ORB2DIABAT_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/qmtool.h"

namespace votca {
namespace xtp {

/// Exports BSE states for the external diabatization code Diabat 2.0:
/// `job_name.fchk` (BSE orbitals, Diabat `mofile`) and
/// `job_name_singlet.diabat` / `job_name_triplet.diabat` (Diabat native
/// `statefile`).
class Orb2Diabat final : public QMTool {
 public:
  Orb2Diabat() = default;
  ~Orb2Diabat() = default;

  std::string Identify() const { return "orb2diabat"; }

 protected:
  void ParseOptions(const tools::Property& options);
  bool Run();

 private:
  std::string orbfile_;
  std::string spin_;
  Index nstates_ = 0;
  Logger log_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ORB2DIABAT_H
