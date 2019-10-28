///*
// *            Copyright 2009-2019 The VOTCA Development Team
// *                       (http://www.votca.org)
// *
// *      Licensed under the Apache License, Version 2.0 (the "License")
// *
// * You may not use this file except in compliance with the License.
// * You may obtain a copy of the License at
// *
// *              http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS,
// * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// */

#pragma once
#ifndef VOTCA_XTP_STERNHEIMERW_H
#define VOTCA_XTP_STERNHEIMERW_H

#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/multishift.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/padeapprox.h>

namespace votca {
namespace xtp {
class Orbitals;
class AOBasis;
class SternheimerW {
 public:
  SternheimerW(Orbitals& orbitals, Logger& log)
      : _orbitals(orbitals), _log(log){};

  void Initialize();

  // Calculates the screened coulomb interaction matrix from the dielectric
  // matrix at frequency w
  std::vector<Eigen::MatrixXcd> ScreenedCoulombOS(std::complex<double> w);
  Eigen::MatrixXcd CalculateScreenedCoulombSC(std::complex<double> w);

  std::vector<Eigen::MatrixXcd> DeltaNOS(std::complex<double> w,
                                         std::string gridtype);

  std::vector<Eigen::MatrixXcd> Polarisability(
      const std::vector<std::complex<double>>& grid_w,
      const std::vector<std::complex<double>>& w, std::string gridtype);

  void initializeMultishift(int size);

  void initializePade(int size);

 private:
  Logger& _log;

  Orbitals& _orbitals;

  PadeApprox _pade;

  Multishift _multishift;

  int _num_occ_lvls;

  int _basis_size;

  Eigen::MatrixXcd _Hamiltonian_Matrix;

  Eigen::MatrixXcd _overlap_Matrix;

  Eigen::MatrixXcd _density_Matrix;

  Eigen::MatrixXd _mo_coefficients;
  Eigen::VectorXd _mo_energies;

  // returns the overlap matrix for all occupied states
  Eigen::MatrixXd OverlapMatrix() const;
  // returns the density matrix for all occupied states
  Eigen::MatrixXd DensityMatrix() const;
  // returns the hamiltonian matrix for all occupied states
  Eigen::MatrixXd Hamiltonian() const;
  // Calculates coulomb matrix
  Eigen::MatrixXcd CoulombMatrix(Eigen::Vector3d gridpoint) const;
  // returns the expansion of the (screened) Coulomb interaction operator
  Eigen::MatrixXcd CalculateDeltaVExpansion(Eigen::MatrixXcd deltaV);
  // sets up the left hand side of the sternheimer equation
  Eigen::MatrixXcd SternheimerLHS(const Eigen::MatrixXcd& hamiltonian,
                                  const Eigen::MatrixXcd& overlap, double eps,
                                  std::complex<double> w, bool pm) const;
  // sets up the right hand side of the sternheimer equation
  Eigen::VectorXcd SternheimerRHS(const Eigen::MatrixXcd& overlap,
                                  const Eigen::MatrixXcd& density,
                                  const Eigen::MatrixXcd& pertubation,
                                  const Eigen::VectorXcd& coeff) const;
  // solves the sternheimer equation via the Biconjugate gradient method
  Eigen::VectorXcd SternheimerSolve(const Eigen::MatrixXcd& LHS,
                                    const Eigen::VectorXcd& RHS);

  std::vector<Eigen::MatrixXcd> DeltaNOneShot(std::vector<std::complex<double>> w,
                                           Eigen::Vector3d r)const;
};
}  // namespace xtp
}  // namespace votca
#endif //VOTCA_VOTCA-XTP_