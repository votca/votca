/*
 *            Copyright 2009-2024 The VOTCA Development Team
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
#include "votca/xtp/dftgradient.h"

namespace votca {
namespace xtp {

Eigen::MatrixXd DFTGradient::NuclearRepulsionDerivative(
    const QMMolecule& mol) {
  Index natoms = mol.size();
  Eigen::MatrixXd deriv = Eigen::MatrixXd::Zero(natoms, 3);

  // dE_nn/dR_A = sum_{B != A} Z_A Z_B * d(1/R_AB)/dR_A
  //            = -sum_{B != A} Z_A Z_B (R_A - R_B) / |R_A - R_B|^3
  //
  // NOTE: this returns the GRADIENT dE/dR, not the force -dE/dR -- the
  // sign convention for the final assembled total gradient (which this
  // feeds into) is a decision for whatever code combines this with the
  // RI-J and XC terms, not fixed here. Keep this consistent when
  // validating against finite differences: a finite difference of the
  // energy directly gives dE/dR, matching this function's output as-is,
  // with no extra sign flip needed.
  for (Index a = 0; a < natoms; ++a) {
    double Za = static_cast<double>(mol[a].getNuccharge());
    const Eigen::Vector3d& Ra = mol[a].getPos();
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    for (Index b = 0; b < natoms; ++b) {
      if (b == a) {
        continue;
      }
      double Zb = static_cast<double>(mol[b].getNuccharge());
      const Eigen::Vector3d& Rb = mol[b].getPos();
      Eigen::Vector3d Rab_vec = Ra - Rb;
      double Rab = Rab_vec.norm();
      sum += Za * Zb * Rab_vec / (Rab * Rab * Rab);
    }
    deriv.row(a) = -sum.transpose();
  }
  return deriv;
}

}  // namespace xtp
}  // namespace votca
