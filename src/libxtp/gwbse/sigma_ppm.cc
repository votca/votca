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

Eigen::VectorXd Sigma_PPM::CalcCorrelationDiag(
    const Eigen::VectorXd& frequencies) const {

  Eigen::VectorXd RPAEnergies = _rpa.getRPAInputEnergies();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
  const int levelsum = _Mmn.nsize();  // total number of bands
  const int gwsize = _Mmn.auxsize();  // size of the GW basis
  const int lumo = _opt.homo + 1;
  const int qpmin_offset = _opt.qpmin - _opt.rpamin;
  // loop over all GW levels
#pragma omp parallel for
  for (int gw_level = 0; gw_level < _qptotal; gw_level++) {
    const double qpmin = frequencies(gw_level);
    double sigma_c = 0.0;
    // loop over all functions in GW basis
    for (int i_gw = 0; i_gw < gwsize; i_gw++) {
      // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc
      // PPM_construct_parameters
      if (_ppm.getPpm_weight()(i_gw) < 1.e-9) {
        continue;
      }

      const Eigen::VectorXd Mmn2 =
          _Mmn[gw_level + qpmin_offset].col(i_gw).cwiseAbs2();
      const double ppm_freq = _ppm.getPpm_freq()(i_gw);
      const double fac = 0.5 * _ppm.getPpm_weight()(i_gw) * ppm_freq;
      Eigen::ArrayXd denom = qpmin - RPAEnergies.array();
      denom.segment(0, lumo) += _ppm.getPpm_freq()(i_gw);
      denom.segment(lumo, levelsum - lumo) -= _ppm.getPpm_freq()(i_gw);
      Stabilize(denom);
      sigma_c += fac * (denom.inverse() * Mmn2.array()).sum();

    }  // GW functions
    result(gw_level) = sigma_c;
  }  // all bands
  return result;
}

void Sigma_PPM::Stabilize(Eigen::ArrayXd& denom) const {
  const double fourpi = 4 * boost::math::constants::pi<double>();
  for (int i = 0; i < denom.size(); ++i) {
    if (std::abs(denom[i]) < 0.25) {
      denom[i] = denom[i] / (0.5 * (1.0 - std::cos(fourpi * denom[i])));
    }
  }
}

Eigen::MatrixXd Sigma_PPM::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);

#pragma omp parallel
  {
    const int lumo = _opt.homo + 1;
    const int levelsum = _Mmn.nsize();  // total number of bands
    const int gwsize = _Mmn.auxsize();  // size of the GW basis
    const Eigen::VectorXd ppm_weight = _ppm.getPpm_weight();
    const Eigen::VectorXd ppm_freqs = _ppm.getPpm_freq();
    const Eigen::VectorXd fac = 0.25 * ppm_weight.cwiseProduct(ppm_freqs);
    const int qpmin_offset = _opt.qpmin - _opt.rpamin;

    const Eigen::VectorXd rpaenergies_thread = _rpa.getRPAInputEnergies();
    ;
#pragma omp for schedule(dynamic)
    for (int gw_level1 = 0; gw_level1 < _qptotal; gw_level1++) {
      const Eigen::MatrixXd& Mmn1 = _Mmn[gw_level1 + qpmin_offset];
      const double qpmin1 = frequencies(gw_level1);
      for (int gw_level2 = gw_level1 + 1; gw_level2 < _qptotal; gw_level2++) {
        const Eigen::MatrixXd& Mmn2 = _Mmn[gw_level2 + qpmin_offset];
        const double qpmin2 = frequencies(gw_level2);
        double sigma_c = 0;
        for (int i_gw = 0; i_gw < gwsize; i_gw++) {
          // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc
          // PPM_construct_parameters
          if (ppm_weight(i_gw) < 1.e-9) {
            continue;
          }

          const Eigen::VectorXd Mmn1xMmn2 =
              Mmn1.col(i_gw).cwiseProduct(Mmn2.col(i_gw));
          Eigen::ArrayXd denom1 = rpaenergies_thread;
          denom1.segment(0, lumo) -= ppm_freqs(i_gw);
          denom1.segment(lumo, levelsum - lumo) += ppm_freqs(i_gw);
          Eigen::ArrayXd denom2 = (qpmin2 - denom1);
          Stabilize(denom2);
          denom1 = (qpmin1 - denom1);
          Stabilize(denom1);
          sigma_c +=
              fac(i_gw) *
              ((denom1.inverse() + denom2.inverse()) * Mmn1xMmn2.array()).sum();
        }
        result(gw_level1, gw_level2) = sigma_c;
        result(gw_level2, gw_level1) = sigma_c;
      }  // GW row
    }    // GW col
  }
  return result;
}

}  // namespace xtp
};  // namespace votca
