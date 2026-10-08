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

// Standard includes
#include <iostream>

#include "votca/xtp/ppm.h"
#include "votca/xtp/rpa_uks.h"

namespace votca {
namespace xtp {

namespace {
template <typename RPAType>
Index ConstructPPMParametersImpl(const RPAType& rpa, Eigen::MatrixXd& ppm_phi,
                                 Eigen::VectorXd& ppm_weight,
                                 Eigen::VectorXd& ppm_freq, double screening_r,
                                 double screening_i) {
  SymmetricEigenSystem es =
      SymmetricEigen(rpa.calculate_epsilon_r(screening_r));
  ppm_phi = std::move(es.vectors);

  ppm_weight = 1 - es.values.array().inverse();

  // epsilon at an imaginary frequency is positive definite
  const Eigen::MatrixXd eps_i = rpa.calculate_epsilon_i(screening_i);
  const Eigen::MatrixXd half = eps_i * ppm_phi;
  const Eigen::MatrixXd epsilon_1_inv = InverseSPD(ppm_phi.transpose() * half);

  ppm_freq.resize(es.values.size());
  Index invalid = 0;
#pragma omp parallel for reduction(+ : invalid)
  for (Index i = 0; i < es.values.size(); i++) {
    if (ppm_weight(i) < 1.e-5) {
      ppm_weight(i) = 0.0;
      ppm_freq(i) = 0.5;
      continue;
    } else {
      double nom = epsilon_1_inv(i, i) - 1.0;
      double frac =
          -1.0 * nom / (nom + ppm_weight(i)) * screening_i * screening_i;
      // frac < 0: this mode has no real plasmon-pole frequency at the
      // imaginary fitting point; |frac| is used, as before, and counted
      if (frac < 0.0) {
        ++invalid;
      }
      ppm_freq(i) = std::sqrt(std::abs(frac));
    }
  }
  return invalid;
}
}  // namespace

void PPM::PPM_construct_parameters(const RPA& rpa) {
  invalid_modes_ = ConstructPPMParametersImpl(
      rpa, ppm_phi_, ppm_weight_, ppm_freq_, screening_r, screening_i);
}

void PPM::PPM_construct_parameters(const RPA_UKS& rpa) {
  invalid_modes_ = ConstructPPMParametersImpl(
      rpa, ppm_phi_, ppm_weight_, ppm_freq_, screening_r, screening_i);
}

}  // namespace xtp
}  // namespace votca