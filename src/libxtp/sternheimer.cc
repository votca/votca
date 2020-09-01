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
 *
 */

#include "votca/xtp/ERIs.h"
#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aomatrix3d.h>
#include <votca/xtp/aopotential.h>
#include <votca/xtp/gauss_hermite_quadrature_constants.h>
#include <votca/xtp/gauss_laguerre_quadrature_constants.h>
#include <votca/xtp/gauss_legendre_quadrature_constants.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/multishift.h>
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

void Sternheimer::initializeMultishift(Index size) {
  _multishift.setMatrixSize(size);
}

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

  double stepsize = (omega_end - omega_start) / steps;
  const double ev2hrt = votca::tools::conv::ev2hrt;
  std::complex<double> d(1, 0);
  std::complex<double> i(0, 1);

  for (Index n = 0; n <= steps; n++) {
    // Converts from input eV to Hartree
    grid.push_back(omega_start * ev2hrt + n * stepsize * ev2hrt +
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

  double alpha = 1000;

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
      if (true) {
        LHS_P = LHS_P + alpha * _density_Matrix.transpose() * _overlap_Matrix;
        LHS_M = LHS_M + alpha * _density_Matrix.transpose() * _overlap_Matrix;
      }
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
    if (perturbationVectoroutput.size() > _opt.max_mixing_history - 1) {
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
      if (perturbationVectorInput.size() > _opt.max_mixing_history - 1) {
        perturbationVectorInput.erase(perturbationVectorInput.begin());
      }
      perturbationVectorInput.push_back(perturbationUsed);
    }
  }

  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << "Warning: Sternheimer cycle not converged" << flush;
  return delta_n_step_one;
}

Eigen::MatrixXcd Sternheimer::DeltaVfromDeltaN(Eigen::MatrixXcd& deltaN) const{

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

  for (Index m = 1; m < Input.size(); m++) {

    c(m - 1) = (DeltaN - Output.at(Output.size() - 1 - m) +
                Input.at(Input.size() - 1 - m))
                   .cwiseProduct(DeltaN)
                   .sum()
                   .real();
    for (Index j = 1; j < Input.size(); j++) {

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

  for (Index n = 1; n < Input.size(); n++) {

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

  for (Index i = 0; i < frequency_evaluation_grid.size(); i++) {
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
  for (Index n = 0; n < frequency_evaluation_grid.size(); n++) {
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

  for (Index n = 0; n < Polar.size(); n++) {
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
void Sternheimer::printIsotropicAverage(
    std::vector<Eigen::Matrix3cd>& polar) const {
  std::vector<std::complex<double>> grid =
      BuildGrid(_opt.start_frequency_grid, _opt.end_frequency_grid,
                _opt.number_output_grid_points, 0);
  XTP_LOG(Log::error, *_pLog)
      << "\n"
      << "#Freq (ev) \t polarizability_isotropic_average" << flush;
  for (Index i = 0; i < polar.size(); i++) {

    XTP_LOG(Log::error, *_pLog)
        << real(grid.at(i)) * votca::tools::conv::hrt2ev << "\t"
        << real((polar.at(i)(2, 2))) + real(polar.at(i)(1, 1)) +
               real(polar.at(i)(0, 0)) / 3
        << flush;
  }
}
std::vector<double> Sternheimer::getIsotropicAverage(
    std::vector<Eigen::Matrix3cd>& polar) const {

  std::vector<double> iA;

  for (Index i = 0; i < polar.size(); i++) {

    iA.push_back(real((polar.at(i)(2, 2))) + real(polar.at(i)(1, 1)) +
                 real(polar.at(i)(0, 0)) / 3);
  }
  return iA;
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

void Sternheimer::printHellmannFeynmanForces(
    std::vector<Eigen::Vector3cd>& EnergyGrad) const {
  QMMolecule mol = _orbitals.QMAtoms();

  XTP_LOG(Log::error, *_pLog) << "\n"
                              << "#Atom_Type "
                              << "Atom_Index "
                              << "Gradient x y z " << flush;

  for (int i = 0; i < EnergyGrad.size(); i++) {
    XTP_LOG(Log::error, *_pLog)
        << mol.at(i).getElement() << " " << i << " " << EnergyGrad[i][0].real()
        << " " << EnergyGrad[i][1].real() << " " << EnergyGrad[i][2].real()
        << flush;
  }
}

std::complex<double> Sternheimer::KoopmanCorrection(Index n,
                                                    double deltaf_n) const {

  std::complex<double> v_nn = std::complex<double>(0.0, 0.0);
  // Setting up Grid for Fxc functional

  Vxc_Grid grid;
  grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                 _dftbasis);
  Vxc_Potential<Vxc_Grid> Vxcpot(grid);

  Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());

  // Build reference density matrix
  Eigen::MatrixXd N_n =
      _mo_coefficients.col(n) * _mo_coefficients.col(n).transpose();
  Eigen::MatrixXd N_ref = (1.0 - deltaf_n) * N_n;
  N_ref += _density_Matrix.real();

  // Build KI potentials
  //(1)
  Eigen::MatrixXcd vhxc_1 = Vxcpot.IntegrateVXC(_density_Matrix).matrix();
  vhxc_1 += _eris.ContractRightIndecesWithMatrix(_density_Matrix);
  //(2)
  Eigen::MatrixXcd vhxc_2 = Vxcpot.IntegrateVXC(N_ref).matrix();
  vhxc_2 += _eris.ContractRightIndecesWithMatrix(N_ref);
  //(3)
  Eigen::MatrixXcd vhxc_3 = Vxcpot.IntegrateVXC(N_ref - N_n).matrix();
  vhxc_3 += _eris.ContractRightIndecesWithMatrix(N_ref - N_n);
  // constant
  std::complex<double> constant = (N_ref).cwiseProduct(vhxc_2).sum();
  constant -= (N_ref - N_n).cwiseProduct(vhxc_3).sum();
  constant -= (N_n).cwiseProduct(vhxc_2).sum();
  // Evaluate expectation value
  return (N_n).cwiseProduct(-vhxc_1 + vhxc_3).sum() + constant;
}

std::complex<double> Sternheimer::KoopmanRelaxationCoeff(
    Index n, double deltaf_n) const {

  std::complex<double> alpha_n = std::complex<double>(0.0, 0.0);
  // Setting up Grid for Fxc functional

  Vxc_Grid grid;
  grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                 _dftbasis);
  Vxc_Potential<Vxc_Grid> Vxcpot(grid);
  Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());

  // Build MO-specific density

  Eigen::MatrixXcd N_n =
      _mo_coefficients.col(n) * _mo_coefficients.col(n).transpose();

  // Build inital perturbation

  Eigen::MatrixXcd FxcInt_init =
      deltaf_n * Fxc(N_n);  // Vxcpot.IntegrateFXC(_density_Matrix, N_n);
  FxcInt_init += _eris.ContractRightIndecesWithMatrix(N_n);

  // Do Sternheimer
  Eigen::MatrixXcd DeltaN = DeltaNSCSternheimer(0.0, FxcInt_init);
  Eigen::MatrixXcd contract = _eris.ContractRightIndecesWithMatrix(DeltaN);
  Eigen::MatrixXcd FxcInt =
      Fxc(DeltaN);  // Vxcpot.IntegrateFXC(_density_Matrix, DeltaN);
  Eigen::MatrixXcd DeltaV = FxcInt_init + contract + FxcInt;
  // Calculate orbital relaxation coeffs
  alpha_n = (N_n).cwiseProduct(DeltaV).sum();
  alpha_n /= (N_n).cwiseProduct(FxcInt_init).sum();
  return alpha_n;
}

void Sternheimer::printKoopmanRelaxationCoeff(std::complex<double> alpha,
                                              Index n) const {
  XTP_LOG(Log::error, *_pLog) << "\n"
                              << "#Orbital "
                              << "alpha " << flush;
  XTP_LOG(Log::error, *_pLog) << n << " " << alpha << flush;
}

void Sternheimer::printKoopman(std::complex<double> alpha,
                               std::complex<double> correction, Index n) const {
  XTP_LOG(Log::error, *_pLog) << "\n"
                              << "#Orbital "
                              << "alpha "
                              << "correction"
                              << " "
                              << "Product" << flush;
  XTP_LOG(Log::error, *_pLog) << n << " " << alpha << " " << correction << " "
                              << alpha * correction << flush;
}

std::vector<Eigen::Vector3cd> Sternheimer::MOEnergyGradient(Index n,
                                                            Index m) const {

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
void Sternheimer::printMOEnergyGradient(
    std::vector<Eigen::Vector3cd>& EnergyGrad, Index n, Index m) const {
  QMMolecule mol = _orbitals.QMAtoms();
  XTP_LOG(Log::error, *_pLog) << "\n"
                              << "#Atom_Type "
                              << "Atom_Index "
                              << "Gradient x y z " << flush;
  for (int i = 0; i < EnergyGrad.size(); i++) {
    XTP_LOG(Log::error, *_pLog)
        << mol.at(i).getElement() << " " << i << " " << EnergyGrad[i][0].real()
        << " " << EnergyGrad[i][1].real() << " " << EnergyGrad[i][2].real()
        << flush;
  }
}

// GW Sternheimer

Eigen::MatrixXcd Sternheimer::GreensFunctionLHS(std::complex<double> w) const {
  return _Hamiltonian_Matrix - w * _overlap_Matrix;
}

Eigen::MatrixXcd Sternheimer::AnalyticGreensfunction(
    std::complex<double> w) const {
  std::complex<double> eta(0., 0.1 * tools::conv::ev2hrt);
  return GreensFunctionLHS(w + eta).colPivHouseholderQr().solve(
      -1 * Eigen::MatrixXcd::Identity(_basis_size, _basis_size));
}
std::complex<double> Sternheimer::Lorentzian(double center,
                                             std::complex<double> freq) const {
  double gamma = 0.1 * tools::conv::ev2hrt;
  //// We avoid on purpose dividing with 1/pi
  // double enumerator = gamma * gamma;
  // double denominator = (std::pow(abs(freq - center), 2) + gamma * gamma);
  // return (enumerator / denominator);
  std::complex<double> j(0., 1.);
  std::complex<double> lorentzian =
      1. / (freq - center - j * gamma) - 1. / (freq - center + j * gamma);
  return lorentzian;
}

Eigen::MatrixXcd Sternheimer::NonAnalyticGreensfunction(
    std::complex<double> freq) const {
  Eigen::MatrixXcd NonAnalyticGreens =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);

  // 2 i pi factor
  std::complex<double> factor(0., 2. * votca::tools::conv::Pi);

  for (Index v = 0; v < _num_occ_lvls; v++) {
    // If the dirac delta is zero (away from its center) avoid doing the product
    NonAnalyticGreens += Lorentzian(_mo_energies[v], freq) *
                         _mo_coefficients.col(v) *
                         _mo_coefficients.col(v).transpose();
  }

  return factor * NonAnalyticGreens;
}

Eigen::MatrixXcd Sternheimer::GreensFunction(
    std::complex<double> frequency) const {

  Eigen::MatrixXcd GreensFunctionMatrix =
      NonAnalyticGreensfunction(frequency) + AnalyticGreensfunction(frequency);
  return GreensFunctionMatrix;
}
Eigen::MatrixXcd Sternheimer::CoulombMatrix(Eigen::Vector3d gridpoint) const {

  AOMultipole aoesp;
  aoesp.FillPotential(_dftbasis, gridpoint);
  return aoesp.Matrix();
}

Eigen::MatrixXcd Sternheimer::ScreenedCoulomb(
    Eigen::Vector3d gridpoint1, std::complex<double> frequency) const {

  XTP_LOG(Log::debug, *_pLog) << "Called screened Coulomb" << flush;

  Eigen::MatrixXcd coulombmatrix = CoulombMatrix(gridpoint1);
  Eigen::MatrixXcd DeltaN = DeltaNSCSternheimer(frequency, coulombmatrix);
  Eigen::MatrixXcd HartreeInt = _eris.ContractRightIndecesWithMatrix(DeltaN);
  Eigen::MatrixXcd FxcInt = DeltaVfromDeltaN(DeltaN);
  Eigen::MatrixXcd DeltaV;

  DeltaV = HartreeInt + FxcInt;  // Watchout! Here we have subtracted
                                 // the bare coulomb potential

  XTP_LOG(Log::debug, *_pLog) << "Screened Coulomb complete" << flush;

  return DeltaV;
}

Eigen::VectorXd Sternheimer::EvaluateBasisAtPosition(
    const AOBasis& dftbasis, const Eigen::Vector3d& pos) const {

  // get value of orbitals at each gridpoint
  Eigen::VectorXd tmat = Eigen::VectorXd::Zero(dftbasis.AOBasisSize());
  for (const AOShell& shell : dftbasis) {
    const double decay = shell.getMinDecay();
    const Eigen::Vector3d& shellpos = shell.getPos();
    Eigen::Vector3d dist = shellpos - pos;
    double distsq = dist.squaredNorm();
    // if contribution is smaller than -ln(1e-10) = 20.27, calc density
    if ((decay * distsq) < -1.0 * std::log(1e-10)) {
      Eigen::VectorBlock<Eigen::VectorXd> tmat_block =
          tmat.segment(shell.getStartIndex(), shell.getNumFunc());
      shell.EvalAOspace(tmat_block, pos);
    }
  }
  return tmat;
}

Eigen::MatrixXcd Sternheimer::SelfEnergy_at_wp(
    std::complex<double> omega, std::complex<double> omega_p) const {
  // This function calculates GW at w and w_p (i.e before integral over
  // frequencies)
  // It perform the spatial grid integration. The final object is a matrix

  XTP_LOG(Log::debug, *_pLog) << "Called print Self-energy at wp" << flush;

  Vxc_Grid _grid;
  _grid.GridSetup("xxcoarse", _orbitals.QMAtoms(), _dftbasis);

  XTP_LOG(Log::debug, *_pLog) << _grid.getGridSize() << flush;

  Eigen::MatrixXcd GF = AnalyticGreensfunction(omega + omega_p);
  Index nthreads = OPENMP::getMaxThreads();
  XTP_LOG(Log::debug, *_pLog) << "Analytic Green's function complete" << flush;
  Eigen::MatrixXcd sigma =
      Eigen::MatrixXcd::Zero(_density_Matrix.rows(), _density_Matrix.cols());
  XTP_LOG(Log::debug, *_pLog) << "Starting spacial integration" << flush;
  for (Index i = 0; i < _grid.getBoxesSize(); ++i) {
    const GridBox& box = _grid[i];
    if (!box.Matrixsize()) {
      continue;
    }

    std::vector<Eigen::MatrixXcd> sigma_thread = std::vector<Eigen::MatrixXcd>(
        nthreads, Eigen::MatrixXcd::Zero(box.Matrixsize(), box.Matrixsize()));

    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

// iterate over gridpoints
#pragma omp parallel for schedule(guided)
    for (Index p = 0; p < box.size(); p++) {
      Eigen::VectorXd ao = box.CalcAOValues(points[p]);
      Eigen::MatrixXd S = ao * ao.transpose();
      const double weight = weights[p];

      // Evaluate bare and screend coulomb potential at point to evaluate the
      // correlation screened Coulomb potential (W_c = W-v). This is evaluated
      // in DeltaNsc

      Eigen::MatrixXcd GW_c = GF * ScreenedCoulomb(points[p], omega_p);

      sigma_thread[OPENMP::getThreadId()] +=
          weight * S * box.ReadFromBigMatrix(GW_c);
    }

    Eigen::MatrixXcd sigma_box = std::accumulate(
        sigma_thread.begin(), sigma_thread.end(),
        Eigen::MatrixXcd::Zero(box.Matrixsize(), box.Matrixsize()).eval());
    box.AddtoBigMatrix(sigma, sigma_box);
  }
  XTP_LOG(Log::debug, *_pLog) << "Self energy at wp complete" << flush;
  return sigma;
}

Eigen::MatrixXcd Sternheimer::SelfEnergy_at_w(
    std::complex<double> omega) const {
  // This function evaluates the frequency integration over w_p, leaving a
  // function of omega Hermite quadrature of order 12: If you want to use
  // different orders for now you have to change the 12 number in getPoints and
  // getAdaptedWeights (now avail. 8,10,12,14,16,18,20,100)
  XTP_LOG(Log::debug, *_pLog) << "Called Self-energy at w" << flush;
  Index order = 8;

  std::complex<double> delta(0., -1e-3);
  double omega_s = 0.0 * tools::conv::ev2hrt;  // Fixed shift on the real axis
  Eigen::MatrixXcd sigma =
      Eigen::MatrixXcd::Zero(_density_Matrix.cols(), _density_Matrix.cols());

  if (_opt.quadrature_scheme == "laguerre") {
    // Laguerre is from 0 to infinity
    Gauss_Laguerre_Quadrature_Constants glqc;
    Eigen::VectorXd _quadpoints = glqc.getPoints(_opt.quadrature_order);
    Eigen::VectorXd _quadadaptedweights =
        glqc.getAdaptedWeights(_opt.quadrature_order);
    for (Index j = 0; j < _opt.quadrature_order; ++j) {
      std::complex<double> gridpoint(0, _quadpoints(j));
      sigma += _quadadaptedweights(j) * SelfEnergy_at_wp(omega, gridpoint);
    }
    sigma *= 2;  // This because later we multiply sigma by -1/2pi
  } else if (_opt.quadrature_scheme == "legendre") {
    // Modified Legendre is from -1 to 1
    Gauss_Legendre_Quadrature_Constants glqc;
    Eigen::VectorXd _quadpoints = glqc.getPoints(_opt.quadrature_order);
    Eigen::VectorXd _quadadaptedweights =
        glqc.getAdaptedWeights(_opt.quadrature_order);
    double x0 = 0.5;
    for (Index j = 0; j < _opt.quadrature_order; ++j) {
      double exponent = (1 + _quadpoints(j)) / (1 - _quadpoints(j));
      double mod_quadpoints = std::pow(x0, exponent);
      double mod_weights = 2 * x0 * _quadadaptedweights(j) /
                           std::pow(1 - _quadadaptedweights(j), 2);
      std::complex<double> gridpoint(0, mod_quadpoints);
      sigma += mod_weights * SelfEnergy_at_wp(omega, gridpoint);
    }
    sigma *= 2;  // This because later we multiply sigma by -1/2pi
  } else if (_opt.quadrature_scheme == "hermite") {
    // Hermite is from -infty to infty
    Gauss_Hermite_Quadrature_Constants glqc;
    Eigen::VectorXd _quadpoints = glqc.getPoints(_opt.quadrature_order);
    Eigen::VectorXd _quadadaptedweights =
        glqc.getAdaptedWeights(_opt.quadrature_order);
    for (Index j = 0; j < _opt.quadrature_order; ++j) {
      std::complex<double> gridpoint(0, _quadpoints(j));
      sigma += _quadadaptedweights(j) * SelfEnergy_at_wp(omega, gridpoint);
    }
  } else {
    XTP_LOG(Log::error, *_pLog)
        << "There no such a thing as the integration scheme you asked" << flush;
  }
  XTP_LOG(Log::debug, *_pLog) << "Self-energy at w complete" << flush;
  return sigma;
}

Eigen::VectorXcd Sternheimer::SelfEnergy_exchange() const {
  Index n_states = _mo_coefficients.cols();
  Eigen::VectorXcd results = Eigen::VectorXcd::Zero(n_states);
  for (Index n = 0; n < n_states; ++n) {
    for (Index v = 0; v < _num_occ_lvls; ++v) {
      Eigen::MatrixXcd P_1 =
          _mo_coefficients.col(v) * _mo_coefficients.col(n).transpose();
      results(n) += (P_1.transpose().cwiseProduct(
                         _eris.ContractRightIndecesWithMatrix(P_1)))
                        .sum();
    }
  }
  return -results;
}

Eigen::VectorXcd Sternheimer::SelfEnergy_diagonal(
    std::complex<double> omega) const {
  XTP_LOG(Log::debug, *_pLog) << "Called SelfEnergy diagonal" << flush;
  Index n_states = _mo_coefficients.cols();
  Eigen::MatrixXcd selfenergy = SelfEnergy_at_w(omega);
  Eigen::VectorXcd results = Eigen::VectorXcd::Zero(n_states);
  for (Index n = 0; n < n_states; ++n) {
    results(n) = (_mo_coefficients.col(n).cwiseProduct(selfenergy *
                                                       _mo_coefficients.col(n)))
                     .sum();
  }
  std::complex<double> prefactor(-1. / (2 * tools::conv::Pi), 0);  // i/(2 eta)
  XTP_LOG(Log::debug, *_pLog) << "Self Energy diagonal complete" << flush;
  return prefactor * results;
}

Eigen::VectorXd Sternheimer::Intercept() const {
  XTP_LOG(Log::debug, *_pLog) << "Called intercept" << flush;
  Index moEs = _mo_energies.size();
  Eigen::VectorXcd Sigma_x = SelfEnergy_exchange();
  Vxc_Grid grid;
  grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                 _dftbasis);
  Vxc_Potential<Vxc_Grid> Vxcpot(grid);
  Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());
  Eigen::MatrixXcd V_xc = Vxcpot.IntegrateVXC(_density_Matrix).matrix();
  Eigen::VectorXcd V_xc_vec = Eigen::VectorXcd::Zero(moEs);
  for (Index n = 0; n < moEs; n++) {
    V_xc_vec(n) =
        (_mo_coefficients.col(n).cwiseProduct(V_xc * _mo_coefficients.col(n)))
            .sum();
  }
  XTP_LOG(Log::debug, *_pLog) << "Intercept complete" << flush;
  return (_orbitals.MOs().eigenvalues() + Sigma_x - V_xc_vec).real();
}

std::complex<double> Sternheimer::SelfEnergy_cohsex(std::complex<double> omega,
                                                    Index n) const {
  XTP_LOG(Log::debug, *_pLog) << "Called SelfEnergy_cohsex" << flush;
  Vxc_Grid _grid;
  _grid.GridSetup("xxcoarse", _orbitals.QMAtoms(), _dftbasis);
  Index nthreads = OPENMP::getMaxThreads();
  std::complex<double> sigma;
  for (Index v = 0; v < _num_occ_lvls; v++) {
    for (Index i = 0; i < _grid.getBoxesSize(); ++i) {
      const GridBox& box = _grid[i];
      if (!box.Matrixsize()) {
        continue;
      }
      std::vector<std::complex<double>> sigma_thread =
          std::vector<std::complex<double>>(nthreads,
                                            std::complex<double>(0., 0.));
      const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
      const std::vector<double>& weights = box.getGridWeights();
// iterate over gridpoints
#pragma omp parallel for schedule(guided)
      for (Index p = 0; p < box.size(); p++) {
        const double weight = weights[p];
        Eigen::VectorXd tmat = EvaluateBasisAtPosition(_dftbasis, points[p]);
        double psi_n = (_mo_coefficients.col(n).cwiseProduct(tmat)).sum();
        double psi_v = (_mo_coefficients.col(v).cwiseProduct(tmat)).sum();
        Eigen::MatrixXcd W_c =
            ScreenedCoulomb(points[p], omega - _mo_energies(v));
        std::complex<double> sc =
            ((_mo_coefficients.col(v))
                 .cwiseProduct(W_c * _mo_coefficients.col(n)))
                .sum();
        sigma_thread[OPENMP::getThreadId()] += weight * psi_n * psi_v * sc;
      }
      std::complex<double> sigma_box =
          std::accumulate(sigma_thread.begin(), sigma_thread.end(),
                          std::complex<double>(0., 0.));
      sigma -= sigma_box;
    }
  }
  XTP_LOG(Log::debug, *_pLog) << "SelfEnergy cohsex complete" << flush;
  return sigma;
}

void Sternheimer::printGW(Index level) const {

  XTP_LOG(Log::debug, *_pLog) << "Called print GW" << flush;

  Index out_points = _opt.number_output_grid_points;
  Index eval_points = _opt.number_of_frequency_grid_points;

  double omega_start = (_opt.start_frequency_grid) * tools::conv::ev2hrt;
  double omega_end = (_opt.end_frequency_grid) * tools::conv::ev2hrt;

  PadeApprox pade1;

  pade1.initialize(eval_points);

  double steps = 0;
  if (eval_points > 1) {
    steps = (omega_end - omega_start) / eval_points;
  }
  Eigen::VectorXd intercept = Intercept();
  for (int j = 0; j < eval_points; ++j) {
    std::complex<double> w(omega_start + j * steps, 0);
    Eigen::VectorXcd sigma_c = SelfEnergy_diagonal(w);
    std::complex<double> sigma_c_sex = SelfEnergy_cohsex(w, level);
    pade1.addPoint(w, sigma_c(level) + sigma_c_sex + intercept(level));
  }
  XTP_LOG(Log::debug, *_pLog) << "Self Energy evaluation succesfull. Writing "
                                 "output using Padé-Approximation."
                              << flush;
  if (out_points > 1) {
    steps = (omega_end - omega_start) / out_points;
  }

  XTP_LOG(Log::error, *_pLog) << "omega"
                              << "\t"
                              << "Real part"
                              << "\t"
                              << "Imag part" << flush;

  for (int j = 0; j < out_points; ++j) {
    double w = omega_start + j * steps;
    XTP_LOG(Log::error, *_pLog)
        << w << "\t" << pade1.evaluatePoint(w).real() << "\t"
        << pade1.evaluatePoint(w).imag() << flush;
  }
  XTP_LOG(Log::debug, *_pLog) << "Print GW complete" << flush;
}

PadeApprox Sternheimer::getGWPade()const{

Index out_points = _opt.number_output_grid_points;
  Index eval_points = _opt.number_of_frequency_grid_points;

  double omega_start = (_opt.start_frequency_grid) * tools::conv::ev2hrt;
  double omega_end = (_opt.end_frequency_grid) * tools::conv::ev2hrt;

  PadeApprox pade;

  pade.initialize(eval_points);

  double steps = 0;
  if (eval_points > 1) {
    steps = (omega_end - omega_start) / eval_points;
  }
  Eigen::VectorXd intercept = Intercept();
  for (int j = 0; j < eval_points; ++j) {
    std::complex<double> w(omega_start + j * steps, 0);
    Eigen::VectorXcd sigma_c = SelfEnergy_diagonal(w);
    std::complex<double> sigma_c_sex = SelfEnergy_cohsex(w, _opt.level);
    pade.addPoint(w, sigma_c(_opt.level) + sigma_c_sex + intercept(_opt.level));
  }
  return pade;

}


}  // namespace xtp

}  // namespace votca