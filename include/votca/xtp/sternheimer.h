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

#include "votca/xtp/ERIs.h"
#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/padeapprox.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/vxc_grid.h>
#include <votca/xtp/vxc_potential.h>

namespace votca {
namespace xtp {

class Sternheimer {
 public:
  Sternheimer();

  Sternheimer(Orbitals& orbitals, Logger* log)
      : _orbitals(orbitals), _pLog(log){};

  // Calculates and saves all matrices needed to perform Sternheimer
  // calculations from DFT
  void setUpMatrices();

  // Options for Sternheimer
  struct options_sternheimer {

    double start_frequency_grid = 0.0;  // in eV
    double end_frequency_grid = 20;     // in eV
    Index number_of_frequency_grid_points = 30;
    double imaginary_shift_pade_approx = 3;  // in eV
    bool do_precalc_fxc = true;        // not recommended for large molecules
    Index number_output_grid_points = 100;
    std::string numerical_Integration_grid_type =
        "coarse";  // xfine fine medium coarse xcoarse
    double perturbation_strength =
        0.1;  // strength of the electric field for polarizability
    Index max_iterations_sc_sternheimer = 100;
    double tolerance_sc_sternheimer = 10E-9;
    double mixing_constant = 0.5;  // 0<mixing_const<1
    Index max_mixing_history = 10;
    Index level = 0; //energy level for GW
    std::string task="gwsternheimer";
  };

  // Edit Options
  void configurate(const options_sternheimer& opt);

  // Bulids the frequency grid for the polarizability calculations
  // Input values in eV
  std::vector<std::complex<double>> BuildGrid(double omega_start,
                                              double omega_end, Index steps,
                                              double imaginary_shift) const;

  // Calculates the Polarizability Tensor for given frequency grid according to
  // Paper https://journals.aps.org/prb/pdf/10.1103/PhysRevB.89.085129
  std::vector<Eigen::Matrix3cd> Polarisability() const;

  std::vector<Eigen::Vector3cd> EnergyGradient() const;

  std::vector<Eigen::Vector3cd> MOEnergyGradient(Index n, Index m) const;

  // Returns Isotropic Average from Polarizability Tensor
  std::vector<double> getIsotropicAverage(
      std::vector<Eigen::Matrix3cd>& polar) const;

 private:
  Logger* _pLog;

  Orbitals& _orbitals;

  AOBasis _dftbasis;
  AOBasis _auxbasis;

  options_sternheimer _opt;

  Index _num_occ_lvls;

  Index _basis_size;

  Eigen::MatrixXcd _Hamiltonian_Matrix;

  Eigen::MatrixXcd _overlap_Matrix;

  Eigen::MatrixXcd _inverse_overlap;

  Eigen::MatrixXd _density_Matrix;

  Eigen::MatrixXd _mo_coefficients;
  Eigen::VectorXd _mo_energies;

  ERIs _eris;
  Eigen::VectorXd _Fxc_presaved;

  // returns the overlap matrix for all occupied states
  Eigen::MatrixXcd OverlapMatrix();
  // returns the density matrix for all occupied states
  Eigen::MatrixXcd DensityMatrix();
  // returns the hamiltonian matrix for all occupied states
  Eigen::MatrixXcd Hamiltonian();
  // Calculates coulomb matrix
  Eigen::MatrixXcd CoulombMatrix();
  
  // sets up the left hand side of the sternheimer equation
  Eigen::MatrixXcd SternheimerLHS(const Eigen::MatrixXcd& hamiltonian,
                                  const Eigen::MatrixXcd& inverse_overlap,
                                  double eps, std::complex<double> omega,
                                  bool pm) const;

  // sets up the right hand side of the sternheimer equation
  Eigen::VectorXcd SternheimerRHS(const Eigen::MatrixXcd& inverse_overlap,
                                  const Eigen::MatrixXcd& density,
                                  const Eigen::MatrixXcd& pertubation,
                                  const Eigen::VectorXd& coeff) const;
  // Calculates the response of the electron density using the self consistent
  // sternheimer method
  Eigen::MatrixXcd DeltaNSCSternheimer(std::complex<double> w,
                            const Eigen::MatrixXcd& pertubation) const;

  //Calculates the Delta V output from SC sternheimer from the Delta N output
  Eigen::MatrixXcd DeltaVfromDeltaN(Eigen::MatrixXcd& deltaN) const;

  // Anderson Mixing with variable history length
  Eigen::MatrixXcd NPAndersonMixing(std::vector<Eigen::MatrixXcd>& Input,
                                    std::vector<Eigen::MatrixXcd>& Output,
                                    double alpha) const;
  Eigen::MatrixXcd Fxc(Eigen::MatrixXcd deltaN) const;

  std::complex<double> Lorentzian(double center, std::complex<double> freq) const;


};
}  // namespace xtp
}  // namespace votca
#endif /* STERNHEIMER_H */