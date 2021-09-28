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

#include "votca/xtp/ERIs.h"
#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aopotential.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/padeapprox.h>
#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/sternheimer.h>
#include <votca/xtp/vxc_grid.h>
#include <votca/xtp/vxc_potential.h>

using std::flush;
namespace votca {
namespace xtp {

// Setting up everything that is needed for a Sternheimer calculation

void Sternheimer::setUpMatrices() {

  XTP_LOG(Log::debug, *_pLog) << "Setting up basis" << flush;

  this->_dftbasis = _orbitals.SetupDftBasis();
  this->_auxbasis = _orbitals.SetupAuxBasis();

  // saving matrices needed from orbitals
  this->_num_occ_lvls = _orbitals.getNumberOfAlphaElectrons();
  this->_basis_size = _orbitals.getBasisSetSize();
  this->_overlap_Matrix = OverlapMatrix();
  this->_density_Matrix = _orbitals.DensityMatrixGroundState();
  this->_mo_coefficients = _orbitals.MOs().eigenvectors();
  this->_mo_energies = _orbitals.MOs().eigenvalues();
  this->_inverse_overlap = _overlap_Matrix.inverse();
  this->_Hamiltonian_Matrix = Hamiltonian();

  ERIs eris;

  XTP_LOG(Log::debug, *_pLog) << "Setting up ERIS" << flush;

  _eris.Initialize(_dftbasis, _auxbasis);
  if (_opt.do_precalc_fxc == true) {
    Vxc_Grid _grid;
    _grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                    _dftbasis);
    Vxc_Potential<Vxc_Grid> Vxcpot(_grid);
    Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());
    this->_Fxc_presaved = Vxcpot.precalcFXC(_density_Matrix);
    XTP_LOG(Log::error, *_pLog) << "Precalculation of Fxc complete" << flush;
  }
}

void Sternheimer::configurate(const options_sternheimer& opt) { _opt = opt; }

Eigen::MatrixXcd Sternheimer::OverlapMatrix() {
  AOOverlap overlap;
  overlap.Fill(_dftbasis);
  return overlap.Matrix().cast<std::complex<double>>();
}

Eigen::MatrixXcd Sternheimer::DensityMatrix() {
  return _orbitals.DensityMatrixGroundState().cast<std::complex<double>>();
}

Eigen::MatrixXcd Sternheimer::Hamiltonian() {

  // rebuilding Hamitonian from eigenvectors and values
  const Eigen::MatrixXd& mo_coefficients = _orbitals.MOs().eigenvectors();
  return (_overlap_Matrix * mo_coefficients *
          _orbitals.MOs().eigenvalues().asDiagonal() *
          mo_coefficients.transpose() * _overlap_Matrix)
      .cast<std::complex<double>>();
}

std::vector<std::complex<double>> Sternheimer::BuildGrid(
    double omega_start, double omega_end, Index steps,
    double imaginary_shift) const {

  std::vector<std::complex<double>> grid;

  double stepsize = (omega_end - omega_start) / (double) steps;
  const double ev2hrt = votca::tools::conv::ev2hrt;
  std::complex<double> d(1, 0);
  std::complex<double> i(0, 1);

  for (Index n = 0; n <= steps; n++) {
    // Converts from input eV to Hartree
    grid.push_back(omega_start * ev2hrt + (double) n * stepsize * ev2hrt +
                   imaginary_shift * i * ev2hrt);
  }
  // Grid has form: stepsize*n+i*imaginary shift
  return grid;
}

Eigen::MatrixXcd Sternheimer::SternheimerLHS(
    const Eigen::MatrixXcd& hamiltonian,
    const Eigen::MatrixXcd& inverse_overlap, double eps,
    std::complex<double> omega, bool pm) const {
  Eigen::MatrixXcd Identity_cmplx =
      Eigen::MatrixXd::Identity(_basis_size, _basis_size)
          .cast<std::complex<double>>();
  std::complex<double> temp = eps + omega;
  if (pm != true) {
    temp = (eps - omega);
  }
  Eigen::MatrixXcd LHS =
      (inverse_overlap * hamiltonian - temp * Identity_cmplx);
  return LHS;
}

Eigen::VectorXcd Sternheimer::SternheimerRHS(
    const Eigen::MatrixXcd& inverse_overlap, const Eigen::MatrixXcd& density,
    const Eigen::MatrixXcd& pertubation, const Eigen::VectorXd& coeff) const {

  Eigen::VectorXcd RHS = -1 * (inverse_overlap - density) * pertubation * coeff;
  return RHS;
}

// The core Self-consistent Sternheimer Algroithm
Eigen::MatrixXcd Sternheimer::DeltaNSCSternheimer(
    std::complex<double> w, const Eigen::MatrixXcd& perturbation) const {

  XTP_LOG(Log::debug, *_pLog) << "Called delta N SC Sternheimer" << flush;

  // Setting up vectors to store old results for Anderson mixing and initial
  // perturbation
  std::vector<Eigen::MatrixXcd> perturbationVectorInput;
  std::vector<Eigen::MatrixXcd> perturbationVectoroutput;

  perturbationVectorInput.push_back(perturbation);
  Eigen::MatrixXcd perturbationUsed = perturbationVectorInput.back();

  Eigen::MatrixXcd delta_n_out_new =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  Eigen::MatrixXcd delta_n_out_old =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  Eigen::MatrixXcd delta_n_step_one =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);

  Vxc_Grid _grid;
  _grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                  _dftbasis);
  Vxc_Potential<Vxc_Grid> Vxcpot(_grid);
  Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());

  for (Index n = 0; n < _opt.max_iterations_sc_sternheimer; n++) {
    Eigen::MatrixXcd FxcInt;
    Eigen::MatrixXcd solution_p =
        Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);
    Eigen::MatrixXcd solution_m =
        Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);

    for (Index v = 0; v < _num_occ_lvls; v++) {

      // Building RHS
      Eigen::MatrixXcd RHS =
          SternheimerRHS(_inverse_overlap, _density_Matrix, perturbationUsed,
                         _mo_coefficients.col(v));

      // Building LHS with +/- omega and solving the system
      Eigen::MatrixXcd LHS_P = SternheimerLHS(
          _Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, true);
      Eigen::MatrixXcd LHS_M = SternheimerLHS(
          _Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, false);
      // LHS_P = LHS_P + alpha * _density_Matrix.transpose() * _overlap_Matrix;
      // LHS_M = LHS_M + alpha * _density_Matrix.transpose() * _overlap_Matrix;
      solution_p.col(v) = LHS_P.colPivHouseholderQr().solve(RHS);
      solution_m.col(v) = LHS_M.colPivHouseholderQr().solve(RHS);
    }

    // Saving previous delta n
    delta_n_out_old = delta_n_out_new;
    // Calculating new delta n
    delta_n_out_new =
        2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
            solution_p.transpose() +
        2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
            solution_m.transpose();

    // Perfomring the to four center Integrals to update delta V
    Eigen::MatrixXcd contract =
        _eris.ContractRightIndecesWithMatrix(delta_n_out_new);
    if (_opt.do_precalc_fxc == true) {
      FxcInt = Fxc(delta_n_out_new);
    } else {
      FxcInt = Vxcpot.IntegrateFXC(_density_Matrix, delta_n_out_new);
    }

    // Check if max mixing history is reached and adding new step to history
    if (perturbationVectoroutput.size() > static_cast<long unsigned int>(_opt.max_mixing_history - 1)) {
      perturbationVectoroutput.erase(perturbationVectoroutput.begin());
    }

    perturbationVectoroutput.push_back((perturbation) + contract + FxcInt);
    // perturbationVectoroutput.push_back((perturbation) + contract);

    double diff =
        (perturbationVectorInput.back() - perturbationVectoroutput.back())
            .squaredNorm();
    if (diff < _opt.tolerance_sc_sternheimer) {
      XTP_LOG(Log::debug, *_pLog) << "SC Sternheimer converged" << flush;
      return delta_n_out_new;
    }
    // Mixing if at least in iteration 2
    if (n == 0) {
      perturbationUsed =
          _opt.mixing_constant * perturbationVectoroutput.back() +
          (1 - _opt.mixing_constant) * perturbationVectorInput.back();
      perturbationVectorInput.push_back(perturbationUsed);
    } else {

      perturbationUsed = (NPAndersonMixing(perturbationVectorInput,
                                           perturbationVectoroutput, 0.5));
      if (perturbationVectorInput.size() > static_cast<long unsigned int>(_opt.max_mixing_history - 1)) {
        perturbationVectorInput.erase(perturbationVectorInput.begin());
      }
      perturbationVectorInput.push_back(perturbationUsed);
    }
  }

  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << "Warning: Sternheimer cycle not converged" << flush;
  return delta_n_step_one;
}

Eigen::MatrixXcd Sternheimer::DeltaVfromDeltaN(Eigen::MatrixXcd& deltaN) const {

  if (_opt.do_precalc_fxc == true) {
    return Fxc(deltaN);
  } else {
    Vxc_Grid _grid;
    _grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                    _dftbasis);
    Vxc_Potential<Vxc_Grid> Vxcpot(_grid);
    Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());
    return Vxcpot.IntegrateFXC(_density_Matrix, deltaN);
  }
}

Eigen::MatrixXcd Sternheimer::NPAndersonMixing(
    std::vector<Eigen::MatrixXcd>& Input, std::vector<Eigen::MatrixXcd>& Output,
    double alpha) const {

  Eigen::MatrixXcd DeltaN = Output.back() - Input.back();

  // Building Linear System for Coefficients
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(Input.size() - 1, Input.size() - 1);
  Eigen::VectorXd c = Eigen::VectorXd::Zero(Input.size() - 1);

  for (Index m = 1; m < Index(Input.size()); m++) {

    c(m - 1) = (DeltaN - Output.at(Output.size() - 1 - m) +
                Input.at(Input.size() - 1 - m))
                   .cwiseProduct(DeltaN)
                   .sum()
                   .real();
    for (Index j = 1; j < Index(Input.size()); j++) {

      A(m - 1, j - 1) =
          (DeltaN - Output.at(Output.size() - 1 - m) +
           Input.at(Input.size() - 1 - m))
              .cwiseProduct(DeltaN - Output.at(Output.size() - 1 - j) +
                            Input.at(Input.size() - 1 - j))
              .sum()
              .real();
    }
  }
  // Solving the System to obtain coefficients
  Eigen::VectorXcd coefficients = A.fullPivHouseholderQr().solve(c);

  // Mixing the Potentials

  Eigen::MatrixXcd OutMixed = Output.back();
  Eigen::MatrixXcd InMixed = Input.back();

  for (Index n = 1; n < Index(Input.size()); n++) {

    OutMixed += coefficients(n - 1) * (Output.at(Output.size() - 1 - n) -
                                       Output.at(Output.size() - 1));
    InMixed += coefficients(n - 1) *
               (Input.at(Input.size() - 1 - n) - Input.at(Input.size() - 1));
  }

  // Returning the linear Mix of Input and Output
  return alpha * OutMixed + (1 - alpha) * InMixed;
}

// Using presave fxc
Eigen::MatrixXcd Sternheimer::Fxc(Eigen::MatrixXcd deltaN) const {

  Index vectorSize = (_basis_size * (_basis_size + 1)) / 2;

  Eigen::MatrixXcd Fxc_sum =
      Eigen::MatrixXcd::Zero(deltaN.cols(), deltaN.cols());

  for (Index i = 0; i < _basis_size; i++) {
    Index sum_i = (i * (i + 1)) / 2;
    for (Index j = i; j < _basis_size; j++) {
      Index index_ij = _basis_size * i - sum_i + j;
      Index index_ij_kl_a =
          vectorSize * index_ij - (index_ij * (index_ij + 1)) / 2;
      for (Index k = 0; k < _basis_size; k++) {
        Index sum_k = (k * (k + 1)) / 2;
        for (Index l = k; l < _basis_size; l++) {
          Index index_kl = _basis_size * k - sum_k + l;
          Index index_ij_kl = index_ij_kl_a + index_kl;
          if (index_ij > index_kl) {
            index_ij_kl = vectorSize * index_kl -
                          (index_kl * (index_kl + 1)) / 2 + index_ij;
          }
          if (l == k) {
            Fxc_sum(i, j) += deltaN(k, l) * _Fxc_presaved(index_ij_kl);
          } else {
            Fxc_sum(i, j) +=
                (deltaN(k, l) + deltaN(l, k)) * _Fxc_presaved(index_ij_kl);
          }
        }
      }
      Fxc_sum(j, i) = Fxc_sum(i, j);
    }
  }
  return Fxc_sum;
}

// Polarizability framework
std::vector<Eigen::Matrix3cd> Sternheimer::Polarisability() const {

  std::vector<std::complex<double>> frequency_evaluation_grid = BuildGrid(
      _opt.start_frequency_grid, _opt.end_frequency_grid,
      _opt.number_of_frequency_grid_points, _opt.imaginary_shift_pade_approx);

  std::vector<std::complex<double>> output_grid =
      BuildGrid(_opt.start_frequency_grid, _opt.end_frequency_grid,
                _opt.number_output_grid_points, 0);

  std::vector<Eigen::Matrix3cd> Polar;
  std::vector<Eigen::Matrix3cd> Polar_pade;

  for (Index i = 0; i < Index(frequency_evaluation_grid.size()); i++) {
    Polar.push_back(Eigen::Matrix3cd::Zero());
  }

  PadeApprox pade_1;
  PadeApprox pade_4;
  PadeApprox pade_6;
  pade_1.initialize(4 * frequency_evaluation_grid.size());
  pade_4.initialize(4 * frequency_evaluation_grid.size());
  pade_6.initialize(4 * frequency_evaluation_grid.size());

  AODipole dipole;
  dipole.Fill(_dftbasis);
#pragma omp parallel for
  for (Index n = 0; n < Index(frequency_evaluation_grid.size()); n++) {
    for (Index i = 0; i < 3; i++) {
      Eigen::MatrixXcd delta_n =
          DeltaNSCSternheimer(frequency_evaluation_grid[n],
                              -_opt.perturbation_strength * dipole.Matrix()[i]);
      for (Index j = i; j < 3; j++) {
        Polar[n](i, j) = -(delta_n.cwiseProduct(dipole.Matrix()[j])).sum();
      }
    }
    for (Index i = 2; i < 3; i++) {
      for (Index j = i + 1; j < 3; j++) {
        Polar[n](j, i) = conj(Polar[n](i, j));
      }
    }
  }

  for (Index n = 0; n < Index(Polar.size()); n++) {
    pade_1.addPoint(frequency_evaluation_grid[n], Polar[n](0, 0));
    pade_1.addPoint(conj(frequency_evaluation_grid[n]), conj(Polar[n](0, 0)));
    pade_1.addPoint(-frequency_evaluation_grid[n], conj(Polar[n](0, 0)));
    pade_1.addPoint(-conj(frequency_evaluation_grid[n]), Polar[n](0, 0));

    pade_4.addPoint(frequency_evaluation_grid[n], Polar[n](1, 1));
    pade_4.addPoint(conj(frequency_evaluation_grid[n]), conj(Polar[n](1, 1)));
    pade_4.addPoint(-frequency_evaluation_grid[n], conj(Polar[n](1, 1)));
    pade_4.addPoint(-conj(frequency_evaluation_grid[n]), Polar[n](1, 1));

    pade_6.addPoint(frequency_evaluation_grid[n], Polar[n](2, 2));
    pade_6.addPoint(conj(frequency_evaluation_grid[n]), conj(Polar[n](2, 2)));
    pade_6.addPoint(-frequency_evaluation_grid[n], conj(Polar[n](2, 2)));
    pade_6.addPoint(-conj(frequency_evaluation_grid[n]), Polar[n](2, 2));
  }

  for (std::complex<double> w : output_grid) {
    Polar_pade.push_back(Eigen::Matrix3cd::Zero());
    Polar_pade[Polar_pade.size() - 1](0, 0) = pade_1.evaluatePoint(w);
    Polar_pade[Polar_pade.size() - 1](1, 1) = pade_4.evaluatePoint(w);
    Polar_pade[Polar_pade.size() - 1](2, 2) = pade_6.evaluatePoint(w);
  }
  return Polar_pade;
}

// Energy Gradient framework
std::vector<Eigen::Vector3cd> Sternheimer::EnergyGradient() const {

  QMMolecule mol = _orbitals.QMAtoms();

  // Setting up Grid for Fxc functional

  Index number_of_atoms = mol.size();

  std::vector<Eigen::Vector3cd> EnergyGrad;

  AO3ddipole ao3dDipole;
  // Loop over Nuclei

  for (int k = 0; k < number_of_atoms; k++) {

    ao3dDipole.setCenter(mol.at(k).getPos());
    ao3dDipole.Fill(_dftbasis);

    double sign = 1.0;

    EnergyGrad.push_back(Eigen::Vector3d::Zero());

    for (int a = 0; a < 3; a++) {

      Vxc_Grid _grid;
      _grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                      _dftbasis);
      Vxc_Potential<Vxc_Grid> Vxcpot(_grid);
      Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());
      Eigen::MatrixXcd DeltaN =
          DeltaNSCSternheimer(0.0, sign * ao3dDipole.Matrix()[a]);
      Eigen::MatrixXcd contract = _eris.ContractRightIndecesWithMatrix(DeltaN);
      Eigen::MatrixXcd FxcInt = Vxcpot.IntegrateFXC(_density_Matrix, DeltaN);
      Eigen::MatrixXcd DeltaV =
          sign * ao3dDipole.Matrix()[a] + contract + FxcInt;
      EnergyGrad[k][a] = _density_Matrix.transpose()
                             .cwiseProduct(DeltaV * mol.at(k).getNuccharge())
                             .sum();
    }
    for (int l = 0; l < number_of_atoms; l++) {
      if (l != k) {
        Eigen::Vector3d distance = (mol.at(k).getPos() - mol.at(l).getPos());
        EnergyGrad[k] += mol.at(k).getNuccharge() * mol.at(l).getNuccharge() *
                         distance / std::pow(distance.norm(), 3);
      }
    }
  }
  return EnergyGrad;
}

std::vector<Eigen::Vector3cd> Sternheimer::MOEnergyGradient(Index n) const {

  QMMolecule mol = _orbitals.QMAtoms();

  // Setting up Grid for Fxc functional

  Index number_of_atoms = mol.size();

  std::vector<Eigen::Vector3cd> EnergyGrad;

  AO3ddipole ao3dDipole;
  // Loop over Nuclei

  for (int k = 0; k < number_of_atoms; k++) {

    ao3dDipole.setCenter(mol.at(k).getPos());
    ao3dDipole.Fill(_dftbasis);

    double sign = 1.0;

    EnergyGrad.push_back(Eigen::Vector3d::Zero());

    for (int a = 0; a < 3; a++) {

      Eigen::MatrixXcd DeltaN =
          DeltaNSCSternheimer(0.0, sign * ao3dDipole.Matrix()[a]);
      Eigen::MatrixXcd contract = _eris.ContractRightIndecesWithMatrix(DeltaN);
      Eigen::MatrixXcd FxcInt =
          Fxc(DeltaN);  // Vxcpot.IntegrateFXC(_density_Matrix, DeltaN);
      Eigen::MatrixXcd DeltaV =
          sign * ao3dDipole.Matrix()[a] + contract + FxcInt;
      EnergyGrad[k][a] = _mo_coefficients.col(n).transpose() *
                         (DeltaV * mol.at(k).getNuccharge()) *
                         _mo_coefficients.col(n);
    }
  }
  return EnergyGrad;
}
}  // namespace xtp

}  // namespace votca