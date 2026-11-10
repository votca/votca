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

// Standard includes
#include <vector>

// Local VOTCA includes
#include "orb2diabat.h"
#include "votca/xtp/diabatwriter.h"
#include "votca/xtp/gaussianwriter.h"
#include "votca/xtp/orbitals.h"

namespace votca {
namespace xtp {

void Orb2Diabat::ParseOptions(const tools::Property& options) {
  orbfile_ = job_name_ + ".orb";
  spin_ = options.get(".spin").as<std::string>();
  nstates_ = options.get(".nstates").as<Index>();
}

bool Orb2Diabat::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, log_) << "Loading data from " << orbfile_ << std::flush;
  orbitals.ReadFromCpt(orbfile_);

  std::vector<QMStateType> types;
  if (spin_ == "singlet" || spin_ == "all") {
    types.push_back(QMStateType(QMStateType::Singlet));
  }
  if (spin_ == "triplet" || spin_ == "all") {
    types.push_back(QMStateType(QMStateType::Triplet));
  }

  // MOs and amplitudes must come from the same orbital set: write the
  // orbitals the BSE was solved in (QSGW-rotated orbitals for QSGW)
  GaussianWriter fchk(log_);
  fchk.UseBSEOrbitals(true);
  fchk.WriteFile(job_name_, orbitals);
  if (orbitals.isQSGW()) {
    XTP_LOG(Log::error, log_)
        << "QSGW run: fchk contains the QSGW-rotated orbitals and QP energies"
        << std::flush;
  }

  DiabatWriter writer;
  for (const QMStateType& type : types) {
    const std::string tag =
        (type == QMStateType::Singlet) ? "singlet" : "triplet";
    const std::string statefile = job_name_ + "_" + tag + ".diabat";
    writer.WriteFile(statefile, orbitals, type, nstates_);
    XTP_LOG(Log::error, log_)
        << "Wrote BSE " << tag << " states to " << statefile
        << "\n  Diabat input:\n    $$statepack\n      mofile = \"" << job_name_
        << ".fchk\"\n      statefile = \"" << statefile
        << "\"\n      state_index = 1..N\n    $$" << std::flush;
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
