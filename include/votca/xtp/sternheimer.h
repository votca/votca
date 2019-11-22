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
 */

#pragma once
#ifndef VOTCA_XTP_STERNHEIMER_H
#define VOTCA_XTP_STERNHEIMER_H

#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {
class Orbitals;
class AOBasis;
class Sternheimer {
 public:
  Sternheimer(Orbitals& orbitals, Logger& log)
      : _orbitals(orbitals), _log(log){};

  // Calculates the dielectric matrix via the non selfconsistent Sternheimer
  // equation for frequency w
  
  void Initialize();
  
  std::vector<Eigen::Matrix3cd> Polarisability(
      double omega_start, double omega_end, int steps, int imaginary_shift, int resolution_output);
  void printIsotropicAverage(std::vector<Eigen::Matrix3cd> polar, std::vector<std::complex<double>> grid);

 private:
  Logger& _log;

  Orbitals& _orbitals;

  Multishift _multishift;

  int _num_occ_lvls;

  int _basis_size;

  Eigen::MatrixXcd _Hamiltonian_Matrix;

  Eigen::MatrixXcd _overlap_Matrix;
  
  Eigen::MatrixXcd _inverse_overlap;

  Eigen::MatrixXcd _density_Matrix;

  Eigen::MatrixXd _mo_coefficients;
  Eigen::VectorXd _mo_energies;

  void initializeMultishift(int size);

  void initializePade(int size);
  // returns the overlap matrix for all occupied states
  Eigen::MatrixXcd OverlapMatrix();
  // returns the density matrix for all occupied states
  Eigen::MatrixXcd DensityMatrix();
  // returns the hamiltonian matrix for all occupied states
  Eigen::MatrixXcd Hamiltonian();
  // Calculates coulomb matrix
  Eigen::MatrixXcd CoulombMatrix();
  
  std::vector<std::complex<double>> BuildGrid(double omega_start, double omega_end, int steps, int imaginary_shift);
  
  
  std::vector<Eigen::MatrixXcd> DipoleIntegral(); 
  // sets up the left hand side of the sternheimer equation
  Eigen::MatrixXcd SternheimerLHS(const Eigen::MatrixXcd& hamiltonian,
                                  const Eigen::MatrixXcd& inverse_overlap,
                                  double eps,
                                  std::complex<double> omega,
                                  bool pm) const;
  // sets up the right hand side of the sternheimer equation
  Eigen::VectorXcd SternheimerRHS(const Eigen::MatrixXcd& inverse_overlap,
                                  const Eigen::MatrixXcd& density,
                                  const Eigen::MatrixXcd& pertubation,
                                  const Eigen::VectorXcd& coeff) const;
  // Calculates the pertubation of the electron density using the non
  // self-consistent one shot approach
  std::vector<Eigen::MatrixXcd> DeltaNOneShot(std::vector<std::complex<double>> w,
                                               const Eigen::MatrixXd& pertubation)const;
  Eigen::MatrixXcd DeltaNSelfConsistent(std::complex<double> w,
                                               const Eigen::MatrixXd& initGuess)const;
  // Pade Appoximation on complex grid
};
}  // namespace xtp
}  // namespace votca
#endif /* STERNHEIMER_H */