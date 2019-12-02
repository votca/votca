/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#include <votca/xtp/ppm.h>
#include <votca/xtp/sigma_ppm.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

void Sigma_PPM::PrepareScreening() {
  _ppm.PPM_construct_parameters(_rpa);
  _Mmn.MultiplyRightWithAuxMatrix(_ppm.getPpm_phi());
}

double Sigma_PPM::CalcCorrelation(Index gw_level, double frequency) const {
  const Index lumo = _opt.homo + 1;
  const Index levelsum = _Mmn.nsize();   // total number of bands
  const Index auxsize = _Mmn.auxsize();  // size of the GW basis
  const Eigen::VectorXd ppm_weight = _ppm.getPpm_weight();
  const Eigen::VectorXd ppm_freqs = _ppm.getPpm_freq();
  const Index qpmin_offset = _opt.qpmin - _opt.rpamin;
  const Eigen::VectorXd RPAEnergies = _rpa.getRPAInputEnergies();
  double sigma_c = 0;
  for (Index i_aux = 0; i_aux < auxsize; i_aux++) {
    // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc
    // PPM_construct_parameters
    if (ppm_weight(i_aux) < 1.e-9) {
      continue;
    }
    const double ppm_freq = ppm_freqs(i_aux);
    const double fac = 0.5 * ppm_weight(i_aux) * ppm_freq;
    const Eigen::VectorXd Mmn2 =
        _Mmn[gw_level + qpmin_offset].col(i_aux).cwiseAbs2();
    Eigen::ArrayXd denom = frequency - RPAEnergies.array();
    denom.segment(0, lumo) += ppm_freq;
    denom.segment(lumo, levelsum - lumo) -= ppm_freq;
    Stabilize(denom);
    sigma_c += fac * (denom.inverse() * Mmn2.array()).sum();
  }
  return sigma_c;
}

double Sigma_PPM::CalcCorrelation(Index gw_level1, Index gw_level2,
                                  double frequency1, double frequency2) const {
  const Index lumo = _opt.homo + 1;
  const Index levelsum = _Mmn.nsize();   // total number of bands
  const Index auxsize = _Mmn.auxsize();  // size of the GW basis
  const Eigen::VectorXd ppm_weight = _ppm.getPpm_weight();
  const Eigen::VectorXd ppm_freqs = _ppm.getPpm_freq();
  const Index qpmin_offset = _opt.qpmin - _opt.rpamin;
  const Eigen::VectorXd RPAEnergies = _rpa.getRPAInputEnergies();
  double sigma_c = 0;
  for (Index i_aux = 0; i_aux < auxsize; i_aux++) {
    // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc
    // PPM_construct_parameters
    if (ppm_weight(i_aux) < 1.e-9) {
      continue;
    }
    const double ppm_freq = ppm_freqs(i_aux);
    const double fac = 0.25 * ppm_weight(i_aux) * ppm_freq;
    const Eigen::MatrixXd& Mmn1 = _Mmn[gw_level1 + qpmin_offset];
    const Eigen::MatrixXd& Mmn2 = _Mmn[gw_level2 + qpmin_offset];
    const Eigen::VectorXd Mmn1xMmn2 =
        Mmn1.col(i_aux).cwiseProduct(Mmn2.col(i_aux));
    Eigen::ArrayXd denom1 = RPAEnergies;
    denom1.segment(0, lumo) -= ppm_freq;
    denom1.segment(lumo, levelsum - lumo) += ppm_freq;
    Eigen::ArrayXd denom2 = (frequency2 - denom1);
    Stabilize(denom2);
    denom1 = (frequency1 - denom1);
    Stabilize(denom1);
    sigma_c +=
        fac * ((denom1.inverse() + denom2.inverse()) * Mmn1xMmn2.array()).sum();
  }
  return sigma_c;
}

void Sigma_PPM::Stabilize(Eigen::ArrayXd& denom) const {
  const double fourpi = 4 * boost::math::constants::pi<double>();
  for (Index i = 0; i < denom.size(); ++i) {
    if (std::abs(denom[i]) < 0.25) {
      denom[i] = denom[i] / (0.5 * (1.0 - std::cos(fourpi * denom[i])));
    }
  }
}

}  // namespace xtp
}  // namespace votca
