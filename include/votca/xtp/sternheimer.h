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
#include <votca/xtp/multishift.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

class Sternheimer {
 public:
  Sternheimer(Orbitals& orbitals, Logger& log)
      : _orbitals(orbitals), _log(log){};

  // Calculates and saves all matrices needed to perform Sternheimer
  // calculations from DFT
  void Initialize();

  // Calculates the Polarizability Tensor for given frequency grid according to
  // Paper https://journals.aps.org/prb/pdf/10.1103/PhysRevB.89.085129
  std::vector<Eigen::Matrix3cd> Polarisability(double omega_start,
                                               double omega_end, int steps,
                                               double imaginary_shift,
                                               double lorentzian_broadening,
                                               int resolution_output) const;
  // Prints the isotropic average of the polarizability tensor
  void printIsotropicAverage(std::vector<Eigen::Matrix3cd>& polar,
                             std::vector<std::complex<double>>& grid) const;

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

  // Sets up the Multishift solver for linear systems of given size
  void initializeMultishift(int size);

  // Sets up the N-Point Pade approximation
  void initializePade(int size);
  // returns the overlap matrix for all occupied states
  Eigen::MatrixXcd OverlapMatrix();
  // returns the density matrix for all occupied states
  Eigen::MatrixXcd DensityMatrix();
  // returns the hamiltonian matrix for all occupied states
  Eigen::MatrixXcd Hamiltonian();
  // Calculates coulomb matrix
  Eigen::MatrixXcd CoulombMatrix();

  // Bulids the frequency grid for the polarizability calculations
  // Input values in eV
  std::vector<std::complex<double>> BuildGrid(double omega_start,
                                              double omega_end, int steps,
                                              double imaginary_shift) const;

  // Computes the Dipole Integral
  std::vector<Eigen::MatrixXcd> DipoleIntegral();
  // sets up the left hand side of the sternheimer equation
  Eigen::MatrixXcd SternheimerLHS(const Eigen::MatrixXcd& hamiltonian,
                                  const Eigen::MatrixXcd& inverse_overlap,
                                  double eps, std::complex<double> omega,
                                  bool pm) const;

  std::vector<Eigen::MatrixXcd> InvertLHS(std::vector<Eigen::MatrixXcd>& LHS);
  // sets up the right hand side of the sternheimer equation
  Eigen::VectorXcd SternheimerRHS(const Eigen::MatrixXcd& inverse_overlap,
                                  const Eigen::MatrixXcd& density,
                                  const Eigen::MatrixXcd& pertubation,
                                  const Eigen::VectorXd& coeff) const;
  // Calculates the response of the electron density using one shot Sternheimer
  std::vector<Eigen::MatrixXcd> DeltaNOneShot(
      std::vector<std::complex<double>> w,
      const Eigen::MatrixXd& pertubation) const;
  // Calculates the response of the electron density using self consistent
  // Sternheimer
  Eigen::MatrixXcd DeltaNSelfConsistent(std::complex<double> w,
                                        const Eigen::MatrixXd& initGuess) const;
};
}  // namespace xtp
}  // namespace votca
#endif /* STERNHEIMER_H */