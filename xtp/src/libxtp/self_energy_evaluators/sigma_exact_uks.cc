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

#include "votca/xtp/rpa_uks.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"

namespace {

}  // namespace

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
  residues_ = std::vector<Eigen::MatrixXd>(qptotal_);
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; gw_level++) {
    const Index qpoffset = opt_.qpmin - opt_.rpamin;
    const Eigen::MatrixXd& Mmn_i = Mmn_[gw_level + qpoffset];

    Eigen::MatrixXd res =
        Eigen::MatrixXd::Zero(rpatotal_, rpa_omegas_.size());
    for (Index s = 0; s < rpa_omegas_.size(); s++) {
      res.col(s) = Mmn_i * screening_modes_[s];
    }
    residues_[gw_level] = std::move(res);
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


}  // namespace xtp
}  // namespace votca