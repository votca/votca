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

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/globals.h>

// Local VOTCA includes
#include "sigma_ppm.h"
#include "votca/xtp/ppm.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

void Sigma_PPM::PrepareScreening() {
  {
    // the epsilon builds inside are timed separately by the RPA
    OptionalTiming t(timings_, "screening: PPM eigensolver, inverse");
    ppm_.PPM_construct_parameters(rpa_);
  }
  OptionalTiming t(timings_, "screening: rotate Mmn (PPM frame)");
  Mmn_.MultiplyRightWithAuxMatrix(ppm_.getPpm_phi());
}

std::string Sigma_PPM::ScreeningSummary() const {
  // Sigma_c sees only modes with non-zero weight (weights below 1e-5 are
  // set to zero in PPM_construct_parameters).
  const Eigen::ArrayXd w = ppm_.getPpm_weight().array();
  return "PPM modes with weight > 0: " + std::to_string((w > 0.0).count()) +
         " of " + std::to_string(w.size()) +
         " (> 1e-3: " + std::to_string((w > 1e-3).count()) +
         ", > 1e-2: " + std::to_string((w > 1e-2).count()) + ")";
}

// sum over modes P and levels n of
//   fac_P |M(n,P)|^2 term(t),  t = w_k - e_n +- Omega_P
// (+ for occupied, - for virtual n), for every frequency w_k. This runs
// inside the parallel QP solver: the array expressions are evaluated lazily
// inside the sums, so there are no temporaries (allocations per mode and
// frequency, millions per iteration, were slow on macOS), but the sums are
// still vectorised.
template <class Term>
void Sigma_PPM::AccumulateDiag(Index gw_level, const double* freqs, Index nfreq,
                               double* out, Term term) const {
  // first virtual level, counted from rpamin like the energies and the
  // rows of Mmn_
  const Index lumo = opt_.homo + 1 - opt_.rpamin;
  const Index levelsum = Mmn_.nsize();  // total number of bands
  const Index qpmin_offset = opt_.qpmin - opt_.rpamin;
  const Eigen::MatrixXd& M = Mmn_[gw_level + qpmin_offset];
  const auto e = rpa_.getRPAInputEnergies().array();
  const Eigen::VectorXd& weight = ppm_.getPpm_weight();
  const Eigen::VectorXd& omega = ppm_.getPpm_freq();
  for (Index i_aux = 0; i_aux < Mmn_.auxsize(); i_aux++) {
    // the ppm_weights smaller 1.e-5 are set to zero in
    // PPM_construct_parameters
    if (weight(i_aux) < 1.e-9) {
      continue;
    }
    const double ppm_freq = omega(i_aux);
    const double fac = 0.5 * weight(i_aux) * ppm_freq;
    // lazily evaluated, vectorised array expressions: no temporaries
    const auto m2 = M.col(i_aux).array().square();
    for (Index k = 0; k < nfreq; ++k) {
      const auto t_occ = (freqs[k] + ppm_freq) - e.head(lumo);
      const auto t_virt = (freqs[k] - ppm_freq) - e.tail(levelsum - lumo);
      out[k] += fac * ((m2.head(lumo) * term(t_occ)).sum() +
                       (m2.tail(levelsum - lumo) * term(t_virt)).sum());
    }
  }
}

double Sigma_PPM::CalcCorrelationDiagElement(Index gw_level,
                                             double frequency) const {
  CountDiagEval();
  const double eta2 = opt_.eta * opt_.eta;
  double sigma = 0.0;
  AccumulateDiag(gw_level, &frequency, 1, &sigma,
                 [eta2](const auto& t) { return t / (t.square() + eta2); });
  return sigma;
}

double Sigma_PPM::CalcCorrelationDiagElementDerivative(Index gw_level,
                                                       double frequency) const {
  const double eta2 = opt_.eta * opt_.eta;
  double dsigma_domega = 0.0;
  AccumulateDiag(gw_level, &frequency, 1, &dsigma_domega,
                 [eta2](const auto& t) {
                   return (eta2 - t.square()) / (t.square() + eta2).square();
                 });
  return dsigma_domega;
}

Eigen::VectorXd Sigma_PPM::CalcCorrelationDiagElements(
    Index gw_level, const Eigen::VectorXd& frequencies) const {
  const Index nfreq = frequencies.size();
  for (Index i = 0; i < nfreq; ++i) {
    CountDiagEval();
  }
  const double eta2 = opt_.eta * opt_.eta;
  // the integrals of the level are read once for all frequencies
  Eigen::VectorXd sigma = Eigen::VectorXd::Zero(nfreq);
  AccumulateDiag(gw_level, frequencies.data(), nfreq, sigma.data(),
                 [eta2](const auto& t) { return t / (t.square() + eta2); });
  return sigma;
}

double Sigma_PPM::CalcCorrelationOffDiagElement(Index gw_level1,
                                                Index gw_level2,
                                                double frequency1,
                                                double frequency2) const {
  // first virtual level, counted from rpamin like the energies and the
  // rows of Mmn_
  const Index lumo = opt_.homo + 1 - opt_.rpamin;
  const double eta2 = opt_.eta * opt_.eta;
  const Index levelsum = Mmn_.nsize();   // total number of bands
  const Index auxsize = Mmn_.auxsize();  // size of the GW basis
  const Eigen::VectorXd ppm_weight = ppm_.getPpm_weight();
  const Eigen::VectorXd ppm_freqs = ppm_.getPpm_freq();
  const Index qpmin_offset = opt_.qpmin - opt_.rpamin;
  const Eigen::VectorXd RPAEnergies = rpa_.getRPAInputEnergies();
  double sigma_c = 0;
  for (Index i_aux = 0; i_aux < auxsize; i_aux++) {
    // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc
    // PPM_construct_parameters
    if (ppm_weight(i_aux) < 1.e-9) {
      continue;
    }
    const double ppm_freq = ppm_freqs(i_aux);
    const double fac = 0.25 * ppm_weight(i_aux) * ppm_freq;
    const Eigen::MatrixXd& Mmn1 = Mmn_[gw_level1 + qpmin_offset];
    const Eigen::MatrixXd& Mmn2 = Mmn_[gw_level2 + qpmin_offset];
    const Eigen::ArrayXd Mmn1xMmn2 =
        Mmn1.col(i_aux).cwiseProduct(Mmn2.col(i_aux));
    Eigen::ArrayXd temp1 = RPAEnergies;
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

Eigen::MatrixXd Sigma_PPM::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  // first virtual level, counted from rpamin like the energies and the
  // rows of Mmn_
  const Index lumo = opt_.homo + 1 - opt_.rpamin;
  const double eta2 = opt_.eta * opt_.eta;
  const Index levelsum = Mmn_.nsize();
  const Index auxsize = Mmn_.auxsize();
  const Index qpmin_offset = opt_.qpmin - opt_.rpamin;
  const Eigen::VectorXd& ppm_weight = ppm_.getPpm_weight();
  const Eigen::VectorXd& ppm_freqs = ppm_.getPpm_freq();
  const Eigen::ArrayXd energies = rpa_.getRPAInputEnergies().array();
  // Sigma_mn = A_mn + A_nm with A_mn = sum_{l,P} X_m(l,P) M_n(l,P),
  //   X_m(l,P) = fac_P M_m(l,P) t/(t^2+eta^2),  t = w_m - e_l -+ Omega_P,
  // the same terms as CalcCorrelationOffDiagElement, summed differently.
  const Eigen::MatrixXd A =
      PairContraction(levelsum, [&](Index m, Eigen::MatrixXd& X) {
        const Eigen::MatrixXd& M = Mmn_[m + qpmin_offset];
        X.resize(levelsum, auxsize);
        const double w = frequencies(m);
        Eigen::ArrayXd t(levelsum);
        for (Index i_aux = 0; i_aux < auxsize; ++i_aux) {
          if (ppm_weight(i_aux) < 1.e-9) {
            X.col(i_aux).setZero();
            continue;
          }
          const double ppm_freq = ppm_freqs(i_aux);
          const double fac = 0.25 * ppm_weight(i_aux) * ppm_freq;
          t = w - energies;
          t.head(lumo) += ppm_freq;
          t.tail(levelsum - lumo) -= ppm_freq;
          X.col(i_aux) =
              (fac * M.col(i_aux).array() * t / (t.abs2() + eta2)).matrix();
        }
      });
  Eigen::MatrixXd result = A + A.transpose();
  result.diagonal().setZero();
  return result;
}

}  // namespace xtp
}  // namespace votca
