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
#include "votca/xtp/rpa.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"

namespace votca {
namespace xtp {

void Sigma_Exact::PrepareScreening() {
  RPA::rpa_eigensolution rpa_solution = rpa_.Diagonalize_H2p();
  rpa_omegas_ = rpa_solution.omega;
  residues_ = std::vector<Eigen::MatrixXd>(qptotal_);
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; gw_level++) {
    residues_[gw_level] = CalcResidues(gw_level, rpa_solution.XpY);
  }
  return;
}

double Sigma_Exact::CalcCorrelationDiagElement(Index gw_level,
                                               double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  double sigma = 0.0;
  for (Index s = 0; s < rpa_omegas_.size(); s++) {
    const double eigenvalue = rpa_omegas_(s);
    const Eigen::ArrayXd res_12 = residues_[gw_level].col(s).cwiseAbs2();
    Eigen::ArrayXd temp = -rpa_.getRPAInputEnergies().array() + frequency;
    temp.segment(0, n_occ) += eigenvalue;
    temp.segment(n_occ, n_unocc) -= eigenvalue;
    const Eigen::ArrayXd denom = temp.abs2() + eta2;
    sigma += (res_12 * temp / denom).sum();
  }
  return 2 * sigma;
}

double Sigma_Exact::CalcCorrelationDiagElementDerivative(
    Index gw_level, double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  double dsigma_domega = 0.0;
  for (Index s = 0; s < rpa_omegas_.size(); s++) {
    const double eigenvalue = rpa_omegas_(s);
    const Eigen::ArrayXd res_12 = residues_[gw_level].col(s).cwiseAbs2();
    Eigen::ArrayXd temp = -rpa_.getRPAInputEnergies().array() + frequency;
    temp.segment(0, n_occ) += eigenvalue;
    temp.segment(n_occ, n_unocc) -= eigenvalue;
    const Eigen::ArrayXd denom = temp.abs2() + eta2;
    dsigma_domega += ((eta2 - temp.abs2()) * res_12 / denom.abs2()).sum();
  }
  return 2 * dsigma_domega;
}

double Sigma_Exact::CalcCorrelationOffDiagElement(Index gw_level1,
                                                  Index gw_level2,
                                                  double frequency1,
                                                  double frequency2) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Index rpasize = rpa_omegas_.size();
  double sigma_c = 0.0;
  for (Index s = 0; s < rpasize; s++) {
    const double eigenvalue = rpa_omegas_(s);
    const Eigen::VectorXd& res1 = residues_[gw_level1].col(s);
    const Eigen::VectorXd& res2 = residues_[gw_level2].col(s);
    const Eigen::VectorXd res_12 = res1.cwiseProduct(res2);
    Eigen::ArrayXd temp1 = -rpa_.getRPAInputEnergies().array();
    temp1.segment(0, n_occ) += eigenvalue;
    temp1.segment(n_occ, n_unocc) -= eigenvalue;
    const Eigen::ArrayXd temp2 = temp1 + frequency2;
    temp1 += frequency1;
    const Eigen::ArrayXd numer1 = res_12.array() * temp1;
    const Eigen::ArrayXd numer2 = res_12.array() * temp2;
    const Eigen::ArrayXd denom1 = temp1.abs2() + eta2;
    const Eigen::ArrayXd denom2 = temp2.abs2() + eta2;
    sigma_c += 0.5 * ((numer1 / denom1) + (numer2 / denom2)).sum();
  }
  // Multiply with factor 2.0 to sum over both (identical) spin states
  return 2.0 * sigma_c;
}

Eigen::MatrixXd Sigma_Exact::CalcResidues(Index gw_level,
                                          const Eigen::MatrixXd& XpY) const {
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Index rpasize = n_occ * n_unocc;
  const Index qpoffset = opt_.qpmin - opt_.rpamin;
  vc2index vc = vc2index(0, 0, n_unocc);
  const Eigen::MatrixXd& Mmn_i = Mmn_[gw_level + qpoffset];
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(rpatotal_, rpasize);
  for (Index v = 0; v < n_occ; v++) {     // Sum over v
    auto Mmn_v = Mmn_[v].middleRows(n_occ, n_unocc);
    auto fc = Mmn_v * Mmn_i.transpose();  // Sum over chi
    auto XpY_v = XpY.middleRows(vc.I(v, 0), n_unocc);
    res += fc.transpose() * XpY_v;        // Sum over c
  }
  return res;
}

}  // namespace xtp
}  // namespace votca
