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

#include <votca/tools/constants.h>
#include <votca/tools/globals.h>

#include "sigma_ppm_uks.h"
#include "votca/xtp/ppm.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

void Sigma_PPM_UKS::PrepareScreening() {
  if (ppm_ == nullptr) {
    throw std::runtime_error(
        "Sigma_PPM_UKS: shared PPM parameters were not set before "
        "PrepareScreening().");
  }
  Mmn_.MultiplyRightWithAuxMatrix(ppm_->getPpm_phi());
}

// As Sigma_PPM::AccumulateDiag: lazily evaluated array expressions, no
// temporaries per mode (this runs inside the parallel QP solver).
template <class Term>
double Sigma_PPM_UKS::AccumulateDiag(Index gw_level, double frequency,
                                     Term term) const {
  // first virtual level, counted from rpamin like the energies and the
  // rows of Mmn_
  const Index lumo = opt_.homo + 1 - opt_.rpamin;
  const Index levelsum = Mmn_.nsize();
  const Index qpmin_offset = opt_.qpmin - opt_.rpamin;
  const Eigen::MatrixXd& M = Mmn_[gw_level + qpmin_offset];
  const auto e = getSpinRPAInputEnergies().array();
  const Eigen::VectorXd& weight = ppm_->getPpm_weight();
  const Eigen::VectorXd& omega = ppm_->getPpm_freq();
  double result = 0.0;
  for (Index i_aux = 0; i_aux < Mmn_.auxsize(); i_aux++) {
    if (weight(i_aux) < 1.e-9) {
      continue;
    }
    const double ppm_freq = omega(i_aux);
    const double fac = 0.5 * weight(i_aux) * ppm_freq;
    const auto m2 = M.col(i_aux).array().square();
    const auto t_occ = (frequency + ppm_freq) - e.head(lumo);
    const auto t_virt = (frequency - ppm_freq) - e.tail(levelsum - lumo);
    result += fac * ((m2.head(lumo) * term(t_occ)).sum() +
                     (m2.tail(levelsum - lumo) * term(t_virt)).sum());
  }
  return result;
}

double Sigma_PPM_UKS::CalcCorrelationDiagElement(Index gw_level,
                                                 double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  return AccumulateDiag(gw_level, frequency, [eta2](const auto& t) {
    return t / (t.square() + eta2);
  });
}

double Sigma_PPM_UKS::CalcCorrelationDiagElementDerivative(
    Index gw_level, double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  return AccumulateDiag(gw_level, frequency, [eta2](const auto& t) {
    return (eta2 - t.square()) / (t.square() + eta2).square();
  });
}

double Sigma_PPM_UKS::CalcCorrelationOffDiagElement(Index gw_level1,
                                                    Index gw_level2,
                                                    double frequency1,
                                                    double frequency2) const {
  // first virtual level, counted from rpamin like the energies and the
  // rows of Mmn_
  const Index lumo = opt_.homo + 1 - opt_.rpamin;
  const double eta2 = opt_.eta * opt_.eta;
  const Index levelsum = Mmn_.nsize();
  const Index auxsize = Mmn_.auxsize();
  const Eigen::VectorXd ppm_weight = ppm_->getPpm_weight();
  const Eigen::VectorXd ppm_freqs = ppm_->getPpm_freq();
  const Index qpmin_offset = opt_.qpmin - opt_.rpamin;
  const Eigen::VectorXd& energies = getSpinRPAInputEnergies();

  double sigma_c = 0.0;
  for (Index i_aux = 0; i_aux < auxsize; i_aux++) {
    if (ppm_weight(i_aux) < 1.e-9) {
      continue;
    }

    const double ppm_freq = ppm_freqs(i_aux);
    const double fac = 0.25 * ppm_weight(i_aux) * ppm_freq;

    const Eigen::MatrixXd& Mmn1 = Mmn_[gw_level1 + qpmin_offset];
    const Eigen::MatrixXd& Mmn2 = Mmn_[gw_level2 + qpmin_offset];

    const Eigen::ArrayXd Mmn1xMmn2 =
        Mmn1.col(i_aux).cwiseProduct(Mmn2.col(i_aux));

    Eigen::ArrayXd temp1 = energies.array();
    temp1.segment(0, lumo) -= ppm_freq;
    temp1.segment(lumo, levelsum - lumo) += ppm_freq;

    Eigen::ArrayXd temp2 = (frequency2 - temp1);
    temp1 = (frequency1 - temp1);

    sigma_c +=
        fac * ((temp1 / (temp1.abs2() + eta2) + temp2 / (temp2.abs2() + eta2)) *
               Mmn1xMmn2)
                  .sum();
  }
  return sigma_c;
}

}  // namespace xtp
}  // namespace votca