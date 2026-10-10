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

// Local VOTCA includes
#include "sigma_exact.h"
#include "votca/xtp/exact_pole_sums.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"

namespace votca {
namespace xtp {

void Sigma_Exact::PrepareScreening() {
  RPA::rpa_eigensolution rpa_solution;
  {
    OptionalTiming t(timings_, "screening: RPA pair eigenproblem");
    rpa_solution = rpa_.Diagonalize_H2p();
  }
  OptionalTiming t(timings_, "screening: residues");
  rpa_omegas_ = rpa_solution.omega;
  // residues of level m: M_m Z with the screening modes
  // Z = sum_v M_v,virt^T (X+Y)_v, once for all levels
  const Eigen::MatrixXd modes = ScreeningModes(rpa_solution.XpY);
  const Index B = Mmn_.nsize();
  const Index S = rpa_omegas_.size();
  residues_.resize(B * S, qptotal_);
  const Index qpoffset = opt_.qpmin - opt_.rpamin;
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; gw_level++) {
    Eigen::Map<Eigen::MatrixXd>(residues_.col(gw_level).data(), B, S)
        .noalias() = Mmn_[gw_level + qpoffset] * modes;
  }
}

ExactPoles Sigma_Exact::Poles() const {
  return ExactPoles{rpa_.getRPAInputEnergies(), opt_.homo + 1 - opt_.rpamin,
                    rpa_omegas_, opt_.eta};
}

double Sigma_Exact::CalcCorrelationDiagElement(Index gw_level,
                                               double frequency) const {
  CountDiagEval();
  // factor 2: both (identical) spin channels
  return 2.0 * exact_pole_sums::Value(Poles(), Residues(gw_level), frequency);
}

Eigen::VectorXd Sigma_Exact::CalcCorrelationDiagElements(
    Index gw_level, const Eigen::VectorXd& frequencies) const {
  for (Index i = 0; i < frequencies.size(); ++i) {
    CountDiagEval();
  }
  return 2.0 *
         exact_pole_sums::Values(Poles(), Residues(gw_level), frequencies);
}

double Sigma_Exact::CalcCorrelationDiagElementDerivative(
    Index gw_level, double frequency) const {
  return 2.0 *
         exact_pole_sums::Derivative(Poles(), Residues(gw_level), frequency);
}

double Sigma_Exact::CalcCorrelationOffDiagElement(Index gw_level1,
                                                  Index gw_level2,
                                                  double frequency1,
                                                  double frequency2) const {
  // 2 (both spins) x 1/2 (symmetrised in the two frequencies)
  const Eigen::Vector2d w(frequency1, frequency2);
  const ExactPoles poles = Poles();
  const Eigen::VectorXd r12 =
      residues_.col(gw_level1).cwiseProduct(residues_.col(gw_level2));
  const Eigen::VectorXd sums =
      exact_pole_sums::ValuesDirectWeighted(poles, r12.data(), w);
  return sums(0) + sums(1);
}

Eigen::MatrixXd Sigma_Exact::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  return exact_pole_sums::OffDiagonal(Poles(), residues_, frequencies);
}

Eigen::MatrixXd Sigma_Exact::ScreeningModes(const Eigen::MatrixXd& XpY) const {
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Index rpasize = n_occ * n_unocc;
  vc2index vc = vc2index(0, 0, n_unocc);
  // rows (v, c) of the virtual rows of the hole slices
  Eigen::MatrixXd M(rpasize, Mmn_.auxsize());
#pragma omp parallel for schedule(dynamic)
  for (Index v = 0; v < n_occ; v++) {
    M.middleRows(vc.I(v, 0), n_unocc) = Mmn_[v].middleRows(n_occ, n_unocc);
  }
  return M.transpose() * XpY;
}

}  // namespace xtp
}  // namespace votca
