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
    /// traceless Cartesian quadrupole in the convention of
    /// StaticSite::CalculateCartesianMultipole (as AOMultipole), carried as
    /// the Gaussian quadrupole density c sum_ij Theta_ij d_i d_j g_k
    Eigen::Matrix3d quadrupole = Eigen::Matrix3d::Zero();
  };

  /// Width (bohr) with which point multipoles are represented as Gaussians
  /// in PointGradient. Narrower Gaussians lose precision in the integrals
  /// (factors beta^3/4 and up); the difference from the point is O(width^2)
  /// and is extrapolated away (widths w and 2w).
  static constexpr double kPointWidth = 4e-3;

  /// Gradient of E = Tr(D Potential(sources)) with respect to the atoms of
  /// basis (natoms x 3, returned) and the source positions (sources x 3,
  /// site_gradient), at fixed moments and widths (three-centre Coulomb
  /// derivative integrals; needs libint2 ERI3 derivatives).
  static Eigen::MatrixXd Gradient(const AOBasis& basis,
                                  const std::vector<Source>& sources,
                                  const Eigen::MatrixXd& D,
                                  Eigen::MatrixXd& site_gradient);

  /// The same for point multipoles (the widths of sources are ignored),
  /// i.e. for Tr(D V) with V of AOMultipole up to the sign: Gradient with
  /// Gaussians of width kPointWidth and 2 kPointWidth, Richardson
  /// extrapolated to zero width (the error of a width w is O(w^2)).
  static Eigen::MatrixXd PointGradient(const AOBasis& basis,
                                       std::vector<Source> sources,
                                       const Eigen::MatrixXd& D,
                                       Eigen::MatrixXd& site_gradient);

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
