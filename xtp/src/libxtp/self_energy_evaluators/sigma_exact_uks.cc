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

#include "sigma_exact_uks.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>

#include "votca/xtp/exact_pole_sums.h"
#include "votca/xtp/rpa_uks.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"

namespace {}  // namespace

namespace votca {
namespace xtp {

void Sigma_Exact_UKS::PrepareScreening() {

  // Build Coulomb-active screening modes in the auxiliary basis from the full
  // unrestricted H2p eigenvectors. This suppresses numerically dark / spin-like
  // modes that can appear in the explicit alpha+beta transition basis but
  // should not contribute to the screened Coulomb interaction W.
  const Eigen::VectorXd* cached_omegas = nullptr;
  const std::vector<Eigen::VectorXd>* cached_modes = nullptr;
  rpa_.GetCachedScreeningModes(cached_omegas, cached_modes);

  rpa_omegas_ = *cached_omegas;
  screening_modes_ = *cached_modes;
  // residues of level m: M_m times the modes, one GEMM per level
  Eigen::MatrixXd modes(Mmn_.auxsize(), rpa_omegas_.size());
  for (Index s = 0; s < rpa_omegas_.size(); s++) {
    modes.col(s) = screening_modes_[std::size_t(s)];
  }
  residues_ = std::vector<Eigen::MatrixXd>(qptotal_);
  const Index qpoffset = opt_.qpmin - opt_.rpamin;
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; gw_level++) {
    residues_[gw_level] = Mmn_[gw_level + qpoffset] * modes;
  }
}

namespace {
ExactPoles UKSPoles(const Eigen::VectorXd& energies, Index n_occ,
                    const Eigen::VectorXd& omegas, double eta) {
  return ExactPoles{energies, n_occ, omegas, eta};
}
}  // namespace

double Sigma_Exact_UKS::CalcCorrelationDiagElement(Index gw_level,
                                                   double frequency) const {
  const ExactPoles poles =
      UKSPoles(getSpinRPAInputEnergies(), opt_.homo + 1 - opt_.rpamin,
               rpa_omegas_, opt_.eta);
  return exact_pole_sums::Value(poles, residues_[gw_level].data(), frequency);
}

double Sigma_Exact_UKS::CalcCorrelationDiagElementDerivative(
    Index gw_level, double frequency) const {
  const ExactPoles poles =
      UKSPoles(getSpinRPAInputEnergies(), opt_.homo + 1 - opt_.rpamin,
               rpa_omegas_, opt_.eta);
  return exact_pole_sums::Derivative(poles, residues_[gw_level].data(),
                                     frequency);
}

double Sigma_Exact_UKS::CalcCorrelationOffDiagElement(Index gw_level1,
                                                      Index gw_level2,
                                                      double frequency1,
                                                      double frequency2) const {
  const ExactPoles poles =
      UKSPoles(getSpinRPAInputEnergies(), opt_.homo + 1 - opt_.rpamin,
               rpa_omegas_, opt_.eta);
  const Eigen::MatrixXd r12 =
      residues_[gw_level1].cwiseProduct(residues_[gw_level2]);
  const Eigen::VectorXd sums = exact_pole_sums::ValuesDirectWeighted(
      poles, r12.data(), Eigen::Vector2d(frequency1, frequency2));
  return 0.5 * (sums(0) + sums(1));
}

}  // namespace xtp
}  // namespace votca