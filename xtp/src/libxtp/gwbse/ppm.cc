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

// Local VOTCA includes
#include "votca/xtp/ppm.h"

namespace votca {
namespace xtp {

void PPM::PPM_construct_parameters(const RPA& rpa) {

  // Solve Eigensystem
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      rpa.calculate_epsilon_r(screening_r));
  ppm_phi_ = es.eigenvectors();

  // store PPM weights from eigenvalues
  ppm_weight_ = 1 - es.eigenvalues().array().inverse();

  // a) phi^t * epsilon(1) * phi e.g. transform epsilon(1) to the same space as
  // epsilon(0)
  Eigen::MatrixXd ortho =
      ppm_phi_.transpose() * rpa.calculate_epsilon_i(screening_i) * ppm_phi_;
  Eigen::MatrixXd epsilon_1_inv = ortho.inverse();
  // determine PPM frequencies
  ppm_freq_.resize(es.eigenvalues().size());
#pragma omp parallel for
  for (Index i = 0; i < es.eigenvalues().size(); i++) {
    if (ppm_weight_(i) < 1.e-5) {
      ppm_weight_(i) = 0.0;
      ppm_freq_(i) = 0.5;  // Hartree
      continue;
    } else {
      double nom = epsilon_1_inv(i, i) - 1.0;
      double frac =
          -1.0 * nom / (nom + ppm_weight_(i)) * screening_i * screening_i;
      ppm_freq_(i) = std::sqrt(std::abs(frac));
    }
  }
  return;
}

}  // namespace xtp
}  // namespace votca
