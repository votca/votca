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

// Third party includes
#include <boost/format.hpp>
#include <utility>

// Local VOTCA includes
#include "votca/xtp/esp2multipole.h"
#include "votca/xtp/espfit.h"
#include "votca/xtp/nbo.h"
#include "votca/xtp/populationanalysis.h"

namespace votca {
namespace xtp {
using std::flush;

void Esp2multipole::Initialize(tools::Property& options) {
  do_svd_ = false;

  use_mulliken_ = false;
  use_CHELPG_ = false;
  use_lowdin_ = false;
  state_ = options.get(".state").as<QMState>();

  method_ = options.get(".method").as<std::string>();

  if (method_ == "mulliken") {
    use_mulliken_ = true;
  } else if (method_ == "loewdin") {
    use_lowdin_ = true;
  } else if (method_ == "CHELPG") {
    use_CHELPG_ = true;
  }

  if (options.exists(".constraints")) {
    std::vector<tools::Property*> prop_region =
        options.Select(".constraints.region");
    Index index = 0;
    for (tools::Property* prop : prop_region) {
      std::string indices = prop->get("indices").as<std::string>();
      QMFragment<double> reg = QMFragment<double>(index, indices);
      index++;
      reg.value() = prop->get("charge").as<double>();
      regionconstraint_.push_back(reg);
      XTP_LOG(Log::error, log_) << "Fit constrained by Region" << flush;
      XTP_LOG(Log::error, log_) << reg;
    }
    std::vector<tools::Property*> prop_pair =
        options.Select(".constraints.pair");
    for (tools::Property* prop : prop_pair) {
      std::vector<Index> pairvec = prop->as<std::vector<Index>>();
      pairconstraint_.emplace_back(pairvec[0], pairvec[1]);
      XTP_LOG(Log::error, log_)
          << "Charges " << pairvec[0] << " " << pairvec[1]
          << " constrained to be equal." << flush;
    }
  }

  gridsize_ = options.get(".gridsize").as<std::string>();

  if (options.exists(".svd")) {
    do_svd_ = true;
    conditionnumber_ = options.get(".svd.conditionnumber").as<double>();
  }

  return;
}

void Esp2multipole::PrintDipoles(const Orbitals& orbitals,
                                 const StaticSegment& seg) const {
  Eigen::Vector3d classical_dip = seg.CalcDipole();

  XTP_LOG(Log::error, log_)
      << "El Dipole from fitted charges [e*bohr]:\n\t\t"
      << boost::format(
             " dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f") %
             classical_dip.x() % classical_dip.y() % classical_dip.z() %
             classical_dip.squaredNorm()
      << flush;
  Eigen::Vector3d qm_dip = orbitals.CalcElDipole(state_);
  XTP_LOG(Log::error, log_)
      << "El Dipole from exact qm density [e*bohr]:\n\t\t"
      << boost::format(
             " dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f") %
             qm_dip.x() % qm_dip.y() % qm_dip.z() % qm_dip.squaredNorm()
      << flush;
}

StaticSegment Esp2multipole::Extractingcharges(const Orbitals& orbitals) const {
  XTP_LOG(Log::error, log_) << "===== Running on " << OPENMP::getMaxThreads()
                            << " threads ===== " << flush;
  StaticSegment result("result", 0);
  if (use_mulliken_) {
    Mulliken mulliken;
    result = mulliken.CalcChargeperAtom(orbitals, state_);
  } else if (use_lowdin_) {
    Lowdin lowdin;
    result = lowdin.CalcChargeperAtom(orbitals, state_);
  } else if (use_CHELPG_) {
    Espfit esp = Espfit(log_);
    if (pairconstraint_.size() > 0) {
      esp.setPairConstraint(pairconstraint_);
    }
    if (regionconstraint_.size() > 0) {
      esp.setRegionConstraint(regionconstraint_);
    }

    if (do_svd_) {
      esp.setUseSVD(conditionnumber_);
    }
    result = esp.Fit2Density(orbitals, state_, gridsize_);
  }

  PrintDipoles(orbitals, result);
  return result;
}

}  // namespace xtp
}  // namespace votca
