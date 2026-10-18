/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

#pragma once
#ifndef VOTCA_XTP_AOGAUSSIANSITES_H
#define VOTCA_XTP_AOGAUSSIANSITES_H

// Standard includes
#include <vector>

// Local VOTCA includes
#include "aobasis.h"
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Coupling of the QM region with Gaussian charges and dipoles (the
 * induced charges and dipoles of charge-response sites).
 *
 * A site k of width s carries the charge density
 *   q g_k(r) + mu . (2 beta) (r - R_k) g_k(r),
 *   g_k(r) = (beta/pi)^3/2 exp(-beta |r - R_k|^2),  beta = 1/s^2,
 * i.e. a normalized Gaussian charge and the Gaussian dipole -mu . grad g_k
 * (dipole moment mu). The integrals are three-centre Coulomb integrals of
 * an AO pair with an s and a p primitive at R_k (libint2), so the QM->site
 * potential and the site->QM Hamiltonian use the same integrals and the
 * coupling is exactly symmetric. Width 0 is not allowed (use AOMultipole
 * for points).
 */
class AOGaussianSites {
 public:
  struct Source {
    Eigen::Vector3d pos;
    double width;
    double charge;
    Eigen::Vector3d dipole;
  };

  /// <mu| phi |nu> with phi the electrostatic potential of the sources
  /// (positive next to a positive charge); the electrons' Hamiltonian gets
  /// minus this matrix.
  static Eigen::MatrixXd Potential(const AOBasis& basis,
                                   const std::vector<Source>& sources);

  /// For the charge density rho(r) = sum_mn D_mn chi_m(r) chi_n(r) (no
  /// sign: electrons are -D): per point k (rows) the potential of rho seen
  /// by a Gaussian of width widths[k] at points[k] (column 0) and its
  /// gradient with respect to points[k] (columns 1-3, = minus the field).
  static Eigen::MatrixXd PotentialAndGradientAt(
      const AOBasis& basis, const Eigen::MatrixXd& D,
      const std::vector<Eigen::Vector3d>& points,
      const std::vector<double>& widths);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOGAUSSIANSITES_H
