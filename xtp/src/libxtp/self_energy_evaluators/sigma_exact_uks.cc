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

#include "votca/xtp/rpa_uks.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"

namespace votca {
namespace xtp {

void Sigma_Exact_UKS::PrepareScreening() {
  RPA_UKS::rpa_eigensolution rpa_solution = rpa_.Diagonalize_H2p();
  rpa_omegas_ = rpa_solution.omega;

  residues_ = std::vector<Eigen::MatrixXd>(qptotal_);
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; gw_level++) {
    residues_[gw_level] = CalcResidues(gw_level, rpa_solution.XpY);
  }
}

double Sigma_Exact_UKS::CalcCorrelationDiagElement(Index gw_level,
                                                   double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Eigen::VectorXd& energies = getSpinRPAInputEnergies();

  double sigma = 0.0;
  for (Index s = 0; s < rpa_omegas_.size(); s++) {
    const double eigenvalue = rpa_omegas_(s);
    const Eigen::ArrayXd res_12 = residues_[gw_level].col(s).cwiseAbs2();

    Eigen::ArrayXd temp = -energies.array() + frequency;
    temp.segment(0, n_occ) += eigenvalue;
    temp.segment(n_occ, n_unocc) -= eigenvalue;

    const Eigen::ArrayXd denom = temp.abs2() + eta2;
    sigma += (res_12 * temp / denom).sum();
  }
  return sigma;
}

double Sigma_Exact_UKS::CalcCorrelationDiagElementDerivative(
    Index gw_level, double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Eigen::VectorXd& energies = getSpinRPAInputEnergies();

  double dsigma_domega = 0.0;
  for (Index s = 0; s < rpa_omegas_.size(); s++) {
    const double eigenvalue = rpa_omegas_(s);
    const Eigen::ArrayXd res_12 = residues_[gw_level].col(s).cwiseAbs2();

    Eigen::ArrayXd temp = -energies.array() + frequency;
    temp.segment(0, n_occ) += eigenvalue;
    temp.segment(n_occ, n_unocc) -= eigenvalue;

    const Eigen::ArrayXd denom = temp.abs2() + eta2;
    dsigma_domega += ((eta2 - temp.abs2()) * res_12 / denom.abs2()).sum();
  }
  return dsigma_domega;
}

double Sigma_Exact_UKS::CalcCorrelationOffDiagElement(Index gw_level1,
                                                      Index gw_level2,
                                                      double frequency1,
                                                      double frequency2) const {
  const double eta2 = opt_.eta * opt_.eta;
  const Index lumo = opt_.homo + 1;
  const Index n_occ = lumo - opt_.rpamin;
  const Index n_unocc = opt_.rpamax - opt_.homo;
  const Eigen::VectorXd& energies = getSpinRPAInputEnergies();

  double sigma_c = 0.0;
  for (Index s = 0; s < rpa_omegas_.size(); s++) {
    const double eigenvalue = rpa_omegas_(s);
    const Eigen::VectorXd& res1 = residues_[gw_level1].col(s);
    const Eigen::VectorXd& res2 = residues_[gw_level2].col(s);
    const Eigen::VectorXd res_12 = res1.cwiseProduct(res2);

    Eigen::ArrayXd temp1 = -energies.array();
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
  return sigma_c;
}

Eigen::MatrixXd Sigma_Exact_UKS::CalcResidues(Index gw_level,
                                              const Eigen::MatrixXd& XpY) const {
  const Index lumo_alpha = rpa_.getHomoAlpha() + 1;
  const Index lumo_beta = rpa_.getHomoBeta() + 1;

  const Index n_occ_alpha = lumo_alpha - opt_.rpamin;
  const Index n_occ_beta = lumo_beta - opt_.rpamin;
  const Index n_unocc_alpha = opt_.rpamax - rpa_.getHomoAlpha();
  const Index n_unocc_beta = opt_.rpamax - rpa_.getHomoBeta();

  const Index size_alpha = n_occ_alpha * n_unocc_alpha;
  const Index size_beta = n_occ_beta * n_unocc_beta;
  const Index rpasize = size_alpha + size_beta;

  const Index qpoffset = opt_.qpmin - opt_.rpamin;
  const Eigen::MatrixXd& Mmn_i = Mmn_[gw_level + qpoffset];

  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(rpatotal_, rpasize);

  vc2index vc_alpha(0, 0, n_unocc_alpha);
  for (Index v = 0; v < n_occ_alpha; v++) {
    auto Mmn_v = Mmn_spin_.alpha[v].middleRows(n_occ_alpha, n_unocc_alpha);
    auto fc = Mmn_v * Mmn_i.transpose();
    auto XpY_v = XpY.middleRows(vc_alpha.I(v, 0), n_unocc_alpha);
    res += fc.transpose() * XpY_v;
  }

  vc2index vc_beta(0, 0, n_unocc_beta);
  for (Index v = 0; v < n_occ_beta; v++) {
    auto Mmn_v = Mmn_spin_.beta[v].middleRows(n_occ_beta, n_unocc_beta);
    auto fc = Mmn_v * Mmn_i.transpose();
    auto XpY_v = XpY.middleRows(size_alpha + vc_beta.I(v, 0), n_unocc_beta);
    res += fc.transpose() * XpY_v;
  }

  return res;
}

}  // namespace xtp
}  // namespace votca