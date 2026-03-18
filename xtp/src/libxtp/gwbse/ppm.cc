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
void ConstructPPMParametersImpl(const RPAType& rpa, Eigen::MatrixXd& ppm_phi,
                                Eigen::VectorXd& ppm_weight,
                                Eigen::VectorXd& ppm_freq, double screening_r,
                                double screening_i) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      rpa.calculate_epsilon_r(screening_r));
  ppm_phi = es.eigenvectors();

  ppm_weight = 1 - es.eigenvalues().array().inverse();

  Eigen::MatrixXd ortho =
      ppm_phi.transpose() * rpa.calculate_epsilon_i(screening_i) * ppm_phi;
  Eigen::MatrixXd epsilon_1_inv = ortho.inverse();

  ppm_freq.resize(es.eigenvalues().size());
#pragma omp parallel for
  for (Index i = 0; i < es.eigenvalues().size(); i++) {
    if (ppm_weight(i) < 1.e-5) {
      ppm_weight(i) = 0.0;
      ppm_freq(i) = 0.5;
      continue;
    } else {
      double nom = epsilon_1_inv(i, i) - 1.0;
      double frac =
          -1.0 * nom / (nom + ppm_weight(i)) * screening_i * screening_i;
      ppm_freq(i) = std::sqrt(std::abs(frac));
    }
  }
}
}  // namespace

void PPM::PPM_construct_parameters(const RPA& rpa) {
  ConstructPPMParametersImpl(rpa, ppm_phi_, ppm_weight_, ppm_freq_, screening_r,
                             screening_i);
}

void PPM::PPM_construct_parameters(const RPA_UKS& rpa) {
  ConstructPPMParametersImpl(rpa, ppm_phi_, ppm_weight_, ppm_freq_, screening_r,
                             screening_i);
}

}  // namespace xtp
}  // namespace votca