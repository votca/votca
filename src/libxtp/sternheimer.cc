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
#include "votca/xtp/adiis.h"
#include "votca/xtp/diis.h"
#include <chrono>
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
namespace votca {
namespace xtp {

void Sternheimer::setUpMatrices() {

  // saving matrices needed from orbitals
  this->_num_occ_lvls = _orbitals.getNumberOfAlphaElectrons();
  this->_basis_size = _orbitals.getBasisSetSize();
  this->_overlap_Matrix = OverlapMatrix();
  this->_density_Matrix = _orbitals.DensityMatrixGroundState();
  this->_mo_coefficients = _orbitals.MOs().eigenvectors();
  this->_mo_energies = _orbitals.MOs().eigenvalues();
  this->_inverse_overlap = _overlap_Matrix.inverse();
  this->_Hamiltonian_Matrix = Hamiltonian();

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();
  ERIs eris;
  _eris.Initialize(dftbasis, auxbasis);

  Vxc_Grid _grid;
  _grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                  dftbasis);
  Vxc_Potential<Vxc_Grid> Vxcpot(_grid);
  Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());
  this->_Fxc_presaved = Vxcpot.precalcFXC(_density_Matrix);
}

void Sternheimer::configurate(const options_sternheimer& opt) { _opt = opt; }

void Sternheimer::initializeMultishift(Index size) {
  _multishift.setMatrixSize(size);
}

Eigen::MatrixXcd Sternheimer::OverlapMatrix() {

  AOBasis basis = _orbitals.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
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

Eigen::MatrixXcd Sternheimer::DeltaNSC(
    std::complex<double> w, const Eigen::MatrixXcd& perturbation) const {

  // auto start = std::chrono::steady_clock::now();

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

  // auto setupinter1 = std::chrono::steady_clock::now();
  // std::cout << "init done: "
  // 	<< std::chrono::duration_cast<std::chrono::milliseconds>(setupinter1 -
  // start).count()
  // 	<< " sec"<<std::endl<<std::endl;

  // Setting up ERIS for four center integral

  // auto setupinter2 = std::chrono::steady_clock::now();
  // std::cout << "ERIS done: "
  // 	<< std::chrono::duration_cast<std::chrono::milliseconds>(setupinter2 -
  // setupinter1).count()
  // 	<< " sec"<<std::endl<<std::endl;

  AOBasis dftbasis = _orbitals.SetupDftBasis();

  // Setting up Grid for Fxc functional
  // Vxc_Grid grid;
  // grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
  //                dftbasis);
  // Vxc_Potential<Vxc_Grid> Vxcpot(grid);
  // Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());

  // auto setupinter3 = std::chrono::steady_clock::now();
  // std::cout << "grid setup done: "
  // 	<< std::chrono::duration_cast<std::chrono::milliseconds>(setupinter3 -
  // setupinter2).count()
  // 	<< " sec"<<std::endl<<std::endl;

  // double alpha = 4*(_mo_energies(_mo_energies.size()-1)-_mo_energies(0));
  double alpha = 1000;
  // Loop until convergence

  // auto inter1 = std::chrono::steady_clock::now();
  // std::cout << "Setup done: "
  // 	<< std::chrono::duration_cast<std::chrono::milliseconds>(inter1 -
  // start).count()
  // 	<< " sec"<<std::endl<<std::endl;

  for (Index n = 0; n < _opt.max_iterations_sc_sternheimer; n++) {
    // auto ref = std::chrono::steady_clock::now();
    // Matrices to store the solutions of the sternheimer equation
    Eigen::MatrixXcd solution_p =
        Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);
    Eigen::MatrixXcd solution_m =
        Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);
    // Loop over all occupied states
    for (Index v = 0; v < _num_occ_lvls; v++) {

      // Building RHS
      Eigen::MatrixXcd RHS =
          SternheimerRHS(_inverse_overlap, _density_Matrix, perturbationUsed,
                         _mo_coefficients.col(v));

      // std::cout<<"RHS "<< RHS <<std::endl;
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

    //   auto inter2 = std::chrono::steady_clock::now();
    // std::cout << "Sternheimer equation solved for all occ state: "
    // 	<< std::chrono::duration_cast<std::chrono::milliseconds>(inter2 -
    // ref).count()
    // 	<< " sec"<<std::endl<<std::endl;

    // Saving previous delta n
    delta_n_out_old = delta_n_out_new;
    // Calculating new delta n
    delta_n_out_new =
        2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
            solution_p.transpose() +
        2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
            solution_m.transpose();

    // auto inter3 = std::chrono::steady_clock::now();
    // std::cout << "Delta N updated: "
    // 	<< std::chrono::duration_cast<std::chrono::milliseconds>(inter3 -
    // inter2).count()
    // 	<< " sec"<<std::endl<<std::endl;
    // Perfomring the to four center Integrals to update delta V
    Eigen::MatrixXcd contract =
        _eris.ContractRightIndecesWithMatrix(delta_n_out_new);

    // auto inter4 = std::chrono::steady_clock::now();
    // std::cout << "Hartree integral done: "
    // 	<< std::chrono::duration_cast<std::chrono::milliseconds>(inter4 -
    // inter3).count()
    // 	<< " sec"<<std::endl<<std::endl;

    // Eigen::MatrixXcd FxcInt =
    //   Vxcpot.IntegrateFXC(_density_Matrix, delta_n_out_new);


    Eigen::MatrixXcd FxcInt = Fxc(delta_n_out_new);


    //auto inter5 = std::chrono::steady_clock::now();
    // std::cout<<"Classic: \n"<<FxcInt<<std::endl<<std::endl;
    // std::cout<<"Presaved: \n"<<FxcInt2<<std::endl<<std::endl;
    // std::cout<<"diff: \n"<<FxcInt-FxcInt2<<std::endl<<std::endl;
    // std::cout<<"diff norm:
    // \n"<<(FxcInt-FxcInt2).norm()<<std::endl<<std::endl;

    // std::cout << "Fxc integral done: "
    // << std::chrono::duration_cast<std::chrono::milliseconds>(inter5 -
    // inter4).count()
    // << " sec"<<std::endl<<std::endl;

    // Check if max mixing history is reached and adding new step to history
    if (perturbationVectoroutput.size() > _opt.max_mixing_history - 1) {
      perturbationVectoroutput.erase(perturbationVectoroutput.begin());
    }

    perturbationVectoroutput.push_back((perturbation) + contract + FxcInt);
    //perturbationVectoroutput.push_back((perturbation) + contract);

    double diff =
        (perturbationVectorInput.back() - perturbationVectoroutput.back())
            .squaredNorm();
    // std::cout << n << " " << diff << std::endl;
    if (diff < _opt.tolerance_sc_sternheimer) {
      //  std::cout << "Frequency: " << w << "Converged after " << n + 1
      //          << " iteration." << std::endl;
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
    // auto inter6 = std::chrono::steady_clock::now();
    // std::cout << "Mixing done, cycle finished: "
    // << std::chrono::duration_cast<std::chrono::milliseconds>(inter6 -
    // inter5).count()
    // << " sec"<<std::endl<<std::endl<<std::endl;
  }

  std::cout << "NOT converged the frequency is w = " << w << std::endl;
  return delta_n_step_one;
}

Eigen::MatrixXcd Sternheimer::AndersonMixing(Eigen::MatrixXcd& inNew,
                                             Eigen::MatrixXcd& inOld,
                                             Eigen::MatrixXcd& outNew,
                                             Eigen::MatrixXcd& outOld,
                                             double alpha) const {

  std::complex<double> beta =
      (outNew - inNew).cwiseProduct(outNew - inNew - outOld + inOld).sum() /
      ((outNew - inNew).cwiseProduct((outOld - inOld))).sum();

  Eigen::MatrixXcd nIn = beta * inOld + (1 - beta) * inNew;
  Eigen::MatrixXcd nOut = beta * outOld + (1 - beta) * outNew;

  return alpha * nOut + (1 - alpha) * nIn;
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

Eigen::MatrixXcd Sternheimer::BroydenMixing(
    std::vector<Eigen::MatrixXcd>& Input, std::vector<Eigen::MatrixXcd>& Output,
    double jacobianScaling) const {

  // Approximation of the first inverse Jacobian
  Eigen::MatrixXcd FirstInverseJacobian =
      jacobianScaling * Eigen::MatrixXcd::Identity(_basis_size, _basis_size);

  Index histSize = Input.size();

  Eigen::MatrixXd gamma = Eigen::MatrixXd::Zero(histSize, histSize);
  Eigen::MatrixXd alpha = Eigen::MatrixXd::Zero(histSize - 1, histSize - 1);
  Eigen::MatrixXd beta = Eigen::MatrixXd::Zero(histSize - 1, histSize - 1);
  Eigen::VectorXd weights = Eigen::VectorXd::Zero(histSize);

  for (Index m = 0; m < histSize; m++) {
    weights(m) = 1 / sqrt(abs((Output.at(m) - Input.at(m))
                                  .cwiseProduct((Output.at(m) - Input.at(m)))
                                  .sum()
                                  .real()));
  }
  for (Index m = 0; m < histSize - 1; m++) {
    for (Index l = 0; l < histSize - 1; l++) {
      alpha(m, l) =
          weights(m) * weights(l) *
          ((Output.at(l + 1) - Input.at(l + 1) - Output.at(l) + Input.at(l)) /
           (Output.at(l + 1) - Input.at(l + 1) - Output.at(l) + Input.at(l))
               .norm())
              .cwiseProduct(((Output.at(m + 1) - Input.at(m + 1) -
                              Output.at(m) + Input.at(m)) /
                             (Output.at(m + 1) - Input.at(m + 1) -
                              Output.at(m) + Input.at(m))
                                 .norm()))
              .sum()
              .real();
    }
  }
  beta = (weights(0) * weights(0) *
              Eigen::MatrixXd::Identity(histSize - 1, histSize - 1) +
          alpha)
             .inverse();
  for (Index m = 0; m < histSize; m++) {
    for (Index l = 0; l < histSize - 1; l++) {
      for (Index k = 0; k < m - 1; k++) {
        gamma(m, l) +=
            weights(k) *
            ((Output.at(k + 1) - Input.at(k + 1) - Output.at(k) + Input.at(k)) /
             (Output.at(k + 1) - Input.at(k + 1) - Output.at(k) + Input.at(k))
                 .norm())
                .cwiseProduct(Input.at(m) - Input.at(m - 1))
                .sum()
                .real() *
            beta(k, l);
      }
    }
  }
  Eigen::MatrixXcd BroydenMix =
      Input.back() + FirstInverseJacobian * (Output.back() - Input.back());
  for (Index n = 0; n < histSize - 1; n++) {

    BroydenMix -=
        weights(n) * gamma(histSize, n) * FirstInverseJacobian *
            ((Output.at(n + 1) - Input.at(n + 1) - Output.at(n) + Input.at(n)) /
             (Output.at(n + 1) - Input.at(n + 1) - Output.at(n) + Input.at(n))
                 .norm()) +
        (Input.at(n + 1) - Input.at(n)) /
            (Output.at(n + 1) - Input.at(n + 1) - Output.at(n) + Input.at(n))
                .norm();
  }

  return BroydenMix;
}

std::vector<Eigen::Matrix3cd> Sternheimer::Polarisability() const {

  std::vector<std::complex<double>> frequency_evaluation_grid = BuildGrid(
      _opt.start_frequency_grid, _opt.end_frequency_grid,
      _opt.number_of_frequency_grid_points, _opt.imaginary_shift_pade_approx);

  std::vector<std::complex<double>> output_grid =
      BuildGrid(_opt.start_frequency_grid, _opt.end_frequency_grid,
                _opt.number_output_grid_points, _opt.lorentzian_broadening);

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

  AOBasis basis = _orbitals.SetupDftBasis();
  AODipole dipole;
  dipole.Fill(basis);
#pragma omp parallel for
  for (Index n = 0; n < frequency_evaluation_grid.size(); n++) {
    for (Index i = 0; i < 3; i++) {
      Eigen::MatrixXcd delta_n =
          DeltaNSC(frequency_evaluation_grid[n],
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
                _opt.number_output_grid_points, _opt.lorentzian_broadening);
  std::cout << "\n"
            << "#Freq (ev) \t polarizability_isotropic_average" << std::endl;
  for (Index i = 0; i < polar.size(); i++) {
    std::cout << real(grid.at(i)) * votca::tools::conv::hrt2ev << "\t"
              << real((polar.at(i)(2, 2))) + real(polar.at(i)(1, 1)) +
                     real(polar.at(i)(0, 0)) / 3
              << std::endl;
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
std::vector<Eigen::Vector3cd> Sternheimer::EnergyGradient() const {

  QMMolecule mol = _orbitals.QMAtoms();

  // Setting up Grid for Fxc functional

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();
  
  Index number_of_atoms = mol.size();

  std::vector<Eigen::Vector3cd> EnergyGrad;

  AO3ddipole ao3dDipole;
  // Loop over Nuclei

  for (int k = 0; k < number_of_atoms; k++) {

    ao3dDipole.setCenter(mol.at(k).getPos());
    ao3dDipole.Fill(dftbasis);

    double sign = 1.0;

    EnergyGrad.push_back(Eigen::Vector3d::Zero());

    for (int a = 0; a < 3; a++) {

      Eigen::MatrixXcd DeltaN = DeltaNSC(0.0, sign * ao3dDipole.Matrix()[a]);
      Eigen::MatrixXcd contract = _eris.ContractRightIndecesWithMatrix(DeltaN);
      Eigen::MatrixXcd FxcInt =
          Fxc(DeltaN);  // Vxcpot.IntegrateFXC(_density_Matrix, DeltaN);
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
  std::cout << "\n"
            << "#Atom_Type "
            << "Atom_Index "
            << "Gradient x y z " << std::endl;
  for (int i = 0; i < EnergyGrad.size(); i++) {
    std::cout << mol.at(i).getElement() << " " << i << " "
              << EnergyGrad[i][0].real() << " " << EnergyGrad[i][1].real()
              << " " << EnergyGrad[i][2].real() << std::endl;
  }
}

std::complex<double> Sternheimer::KoopmanCorrection(Index n,
                                                    double deltaf_n) const {

  std::complex<double> v_nn = std::complex<double>(0.0, 0.0);
  // Setting up Grid for Fxc functional

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();
  // ERIs eris;
  // eris.Initialize(dftbasis, auxbasis);

  Vxc_Grid grid;
  grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                 dftbasis);
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

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();
  // ERIs eris;
  // eris.Initialize(dftbasis, auxbasis);

  Vxc_Grid grid;
  grid.GridSetup(_opt.numerical_Integration_grid_type, _orbitals.QMAtoms(),
                 dftbasis);
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
  Eigen::MatrixXcd DeltaN = DeltaNSC(0.0, FxcInt_init);
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
  std::cout << "\n"
            << "#Orbital "
            << "alpha " << std::endl;
  std::cout << n << " " << alpha << std::endl;
}

void Sternheimer::printKoopman(std::complex<double> alpha,
                               std::complex<double> correction, Index n) const {
  std::cout << "\n"
            << "#Orbital "
            << "alpha "
            << "correction"
            << " "
            << "Product" << std::endl;
  std::cout << n << " " << alpha << " " << correction << " "
            << alpha * correction << std::endl;
}

std::vector<Eigen::Vector3cd> Sternheimer::MOEnergyGradient(Index n,
                                                            Index m) const {

  QMMolecule mol = _orbitals.QMAtoms();

  // Setting up Grid for Fxc functional

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();
  
  Index number_of_atoms = mol.size();

  std::vector<Eigen::Vector3cd> EnergyGrad;

  AO3ddipole ao3dDipole;
  // Loop over Nuclei

  for (int k = 0; k < number_of_atoms; k++) {

    ao3dDipole.setCenter(mol.at(k).getPos());
    ao3dDipole.Fill(dftbasis);

    double sign = 1.0;

    EnergyGrad.push_back(Eigen::Vector3d::Zero());

    for (int a = 0; a < 3; a++) {

      Eigen::MatrixXcd DeltaN = DeltaNSC(0.0, sign * ao3dDipole.Matrix()[a]);
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
  std::cout << "\n"
            << "#Atom_Type "
            << "Atom_Index "
            << "Gradient x y z " << std::endl;
  for (int i = 0; i < EnergyGrad.size(); i++) {
    std::cout << mol.at(i).getElement() << " " << i << " "
              << EnergyGrad[i][0].real() << " " << EnergyGrad[i][1].real()
              << " " << EnergyGrad[i][2].real() << std::endl;
  }
}

Eigen::MatrixXcd Sternheimer::GreensFunctionLHS(std::complex<double> w) const {
  return _Hamiltonian_Matrix - w * _overlap_Matrix;
}

Eigen::MatrixXcd Sternheimer::AnalyticGreensfunction(
    std::complex<double> w) const {
  std::complex<double> eta(0., 0.3 * tools::conv::ev2hrt);
  return GreensFunctionLHS(w + eta).colPivHouseholderQr().solve(
      -1 * Eigen::MatrixXcd::Identity(_basis_size, _basis_size));
}
double Sternheimer::Lorentzian(double center, std::complex<double> freq) const {
  double gamma = 1e-3;
  // We avoid on purpose dividing with 1/pi
  double enumerator = gamma * gamma;
  double denominator = (std::pow(abs(freq - center), 2) + gamma * gamma);
  return (enumerator / denominator);
}

Eigen::MatrixXcd Sternheimer::NonAnalyticGreensfunction(
    std::complex<double> freq) const {
  Eigen::MatrixXcd NonAnalyticGreens =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);

  // 2 i pi factor
  std::complex<double> factor(0., 2. * votca::tools::conv::Pi);

  for (Index v = 0; v < _num_occ_lvls; v++) {
    // If the dirac delta is zero (away from its center) avoid doing the product
    NonAnalyticGreens +=  Lorentzian(_mo_energies[v], freq) *
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

  AOBasis basis = _orbitals.SetupDftBasis();
  AOMultipole aoesp;
  aoesp.FillPotential(basis, gridpoint);
  return aoesp.Matrix();
}

Eigen::MatrixXcd Sternheimer::ScreenedCoulomb(
    Eigen::Vector3d gridpoint1, std::complex<double> frequency) const {

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();

  // Vxc_Potential<Vxc_Grid> Vxcpot(_grid);
  // Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());
  Eigen::MatrixXcd coulombmatrix = CoulombMatrix(gridpoint1);
  Eigen::MatrixXcd DeltaN = DeltaNSC(frequency, coulombmatrix);
  Eigen::MatrixXcd HartreeInt = _eris.ContractRightIndecesWithMatrix(DeltaN);
  Eigen::MatrixXcd FxcInt =
      Fxc(DeltaN);  // Vxcpot.IntegrateFXC(_density_Matrix, DeltaN);
  // Eigen::MatrixXcd DeltaV = coulombmatrix + HartreeInt + FxcInt;
  Eigen::MatrixXcd DeltaV =
      HartreeInt + FxcInt;  // Watchout! Here we have subtracted
                            // the bare coulomb potential
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

Eigen::MatrixXcd Sternheimer::SelfEnergy_at_wp(std::complex<double> omega,
                                               std::complex<double> omega_p) const {
  // This function calculates GW at w and w_p (i.e before integral over
  // frequencies)
  // It perform the spatial grid integration. The final object is a matrix

  AOBasis basis = _orbitals.SetupDftBasis();
  Vxc_Grid _grid;
  _grid.GridSetup("xxcoarse", _orbitals.QMAtoms(), basis);
  //std::cout<<_grid.getGridSize()<<std::endl;

  // std::vector<const Eigen::Vector3d*> gridp = _grid.getGridpoints();

  // for(Index i=0;i<_grid.getGridSize();i++){
  //   Eigen::Vector3d p=*gridp.at(i);
  //   std::cout<<p.transpose()<<std::endl;
  // }

  Eigen::MatrixXcd GF = GreensFunction(omega + omega_p);
  Index nthreads = OPENMP::getMaxThreads();

  Eigen::MatrixXcd sigma =
      Eigen::MatrixXcd::Zero(_density_Matrix.rows(), _density_Matrix.cols());
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

  return sigma;
}

Eigen::MatrixXcd Sternheimer::SelfEnergy_at_wp_regulargrid(
    std::complex<double> omega, double omega_p) const {
  // This function calculates GW at w and w_p (i.e before integral over
  // frequencies)
  // It perform the spatial uniform grid integration. The final object is a
  // matrix
  AOBasis basis = _orbitals.SetupDftBasis();
  Eigen::MatrixXcd GF = GreensFunction(omega + omega_p);
  Index nthreads = OPENMP::getMaxThreads();

  //double _padding = 6.512752;
  double _padding = 0.1;
  Index _xsteps = _opt.gws_grid_spacing;
  Index _ysteps = _opt.gws_grid_spacing;
  Index _zsteps = _opt.gws_grid_spacing;

  const QMMolecule& atoms = _orbitals.QMAtoms();
  std::pair<Eigen::Vector3d, Eigen::Vector3d> minmax =
      atoms.CalcSpatialMinMax();
  double xstart = minmax.first.x() - _padding;
  double xstop = minmax.second.x() + _padding;
  double ystart = minmax.first.y() - _padding;
  double ystop = minmax.second.y() + _padding;
  double zstart = minmax.first.z() - _padding;
  double zstop = minmax.second.z() + _padding;

  double xincr = (xstop - xstart) / double(_xsteps);
  double yincr = (ystop - ystart) / double(_ysteps);
  double zincr = (zstop - zstart) / double(_zsteps);

  double weight = xincr * yincr * zincr;
  double cx = 1;
  double cy = 1;
  double cz = 1;
  std::vector<Eigen::MatrixXcd> sigma_thread = std::vector<Eigen::MatrixXcd>(
      nthreads,
      Eigen::MatrixXcd::Zero(_density_Matrix.rows(), _density_Matrix.cols()));
// eval density at cube grid points
#pragma omp parallel for schedule(guided)
  for (int ix = 0; ix <= _xsteps; ix++) {
    double x = xstart + double(ix) * xincr;
    if (ix == 0) {
      cx = 0.5;
    } else if (ix == _xsteps) {
      cx = 0.5;
    } else {
      cx = 1;
    }
    for (int iy = 0; iy <= _ysteps; iy++) {
      double y = ystart + double(iy) * yincr;
      if (iy == 0) {
        cy = 0.5;
      } else if (iy == _ysteps) {
        cy = 0.5;
      } else {
        cy = 1;
      }

      for (int iz = 0; iz <= _zsteps; iz++) {
        double z = zstart + double(iz) * zincr;
        if (iz == 0) {
          cz = 0.5;
        } else if (iz == _zsteps) {
          cz = 0.5;
        } else {
          cz = 1;
        }
        Eigen::Vector3d pos(x, y, z);
        Eigen::VectorXd tmat = EvaluateBasisAtPosition(basis, pos);
        // Evaluate bare and screend coulomb potential at point to evaluate the
        // correlation screened Coulomb potential (W_c = W-v). This is evaluated
        // in DeltaNsc
        Eigen::MatrixXcd GW_c = GF * ScreenedCoulomb(pos, omega_p);
        sigma_thread[OPENMP::getThreadId()] +=
            cx * cy * cz * weight * tmat * tmat.transpose() * GW_c;
      }  // z-component
    }    // y-component
  }      // x-component

  Eigen::MatrixXcd sigma = std::accumulate(
      sigma_thread.begin(), sigma_thread.end(),
      Eigen::MatrixXcd::Zero(_density_Matrix.rows(), _density_Matrix.cols())
          .eval());
  return sigma;
}

Eigen::MatrixXcd Sternheimer::SelfEnergy_at_w(std::complex<double> omega) const {
  // This function evaluates the frequency integration over w_p, leaving a
  // function of omega Hermite quadrature of order 12: If you want to use
  // different orders for now you have to change the 12 number in getPoints and
  // getAdaptedWeights (now avail. 8,10,12,14,16,18,20,100)
  Index order = 8;

  std::complex<double> delta(0., -1e-3);
  double omega_s = 0.3*tools::conv::ev2hrt; //Fixed shift respect to real axis
  Eigen::MatrixXcd sigma =
      Eigen::MatrixXcd::Zero(_density_Matrix.cols(), _density_Matrix.cols());

  if (_opt.quadrature_scheme == "laguerre") {
    // Laguerre is from 0 to infinity
    Gauss_Laguerre_Quadrature_Constants glqc;
    Eigen::VectorXd _quadpoints = glqc.getPoints(order);
    Eigen::VectorXd _quadadaptedweights = glqc.getAdaptedWeights(order);
    for (Index j = 0; j < order; ++j) {
    std::complex<double> gridpoint(_quadpoints(j),0); 
    sigma += _quadadaptedweights(j) *
             SelfEnergy_at_wp(omega, gridpoint);
  }
  sigma *=2; //This because later we multiply sigma by -1/2pi
  } else if (_opt.quadrature_scheme == "legendre") {
    //Modified Legendre is from -1 to 1
    Gauss_Legendre_Quadrature_Constants glqc;
    Eigen::VectorXd _quadpoints = glqc.getPoints(_opt.quadrature_order);
    Eigen::VectorXd _quadadaptedweights = glqc.getAdaptedWeights(_opt.quadrature_order);
    double x0 = 0.5;
    for (Index j = 0; j < _opt.quadrature_order; ++j) {
      double exponent = (1+_quadpoints(j))/(1-_quadpoints(j));
      double mod_quadpoints = std::pow(x0,exponent);
      double mod_weights =  2*x0*_quadadaptedweights(j)/std::pow(1-_quadadaptedweights(j),2);
    std::complex<double> gridpoint(mod_quadpoints,0); 
    sigma += mod_weights * SelfEnergy_at_wp(omega, gridpoint);
  }
  sigma *=2; //This because later we multiply sigma by -1/2pi
  } else if (_opt.quadrature_scheme == "hermite") {
    //Hermite is from -infty to infty
    Gauss_Hermite_Quadrature_Constants glqc;
    Eigen::VectorXd _quadpoints = glqc.getPoints(_opt.quadrature_order);
    Eigen::VectorXd _quadadaptedweights = glqc.getAdaptedWeights(_opt.quadrature_order);
    for (Index j = 0; j < _opt.quadrature_order; ++j) {
    std::complex<double> gridpoint(_quadpoints(j),0); 
    sigma += _quadadaptedweights(j) *
             SelfEnergy_at_wp(omega, gridpoint);
  }
  } else {
    std::cout << "There no such a thing as the integration scheme you asked"
              << std::endl;
  }
  
  return sigma;
}

Eigen::MatrixXcd Sternheimer::SelfEnergy_at_w_rect(std::complex<double> omega) const {
  // This function evaluates the frequency integration over w_p, using the rectangle rule

  std::vector<double> _quadpoints;
  std::vector<double> _quadadaptedweights;

  for(Index j = 0; j < 25; ++j){

    _quadpoints.push_back(-6.0+j*0.5);
    _quadadaptedweights.push_back(0.5);

  }

  std::complex<double> delta(0., -1e-3);
  Eigen::MatrixXcd sigma =
      Eigen::MatrixXcd::Zero(_density_Matrix.cols(), _density_Matrix.cols());
  for (Index j = 0; j < 25; ++j) {

    Eigen::MatrixXcd SE=SelfEnergy_at_wp(omega, _quadpoints.at(j));

    sigma += _quadadaptedweights.at(j) * SE;

    std::cout<<_quadpoints.at(j)<<" "<<SE.norm()<<std::endl;

  }

  return sigma;
}

Eigen::VectorXcd Sternheimer::SelfEnergy_diagonal(std::complex<double> omega) const {
  Index n_states = _mo_coefficients.cols();
  Eigen::MatrixXcd selfenergy = SelfEnergy_at_w(omega);
  Eigen::VectorXcd results = Eigen::VectorXcd::Zero(n_states);
  for (Index n = 0; n < n_states; ++n) {
    results(n) = (_mo_coefficients.col(n).cwiseProduct(selfenergy *
                                                       _mo_coefficients.col(n)))
                     .sum();
  }
  std::complex<double> prefactor(-1./(2 * tools::conv::Pi),0.);  // i/(2 eta)
  return prefactor * results;
}

Eigen::MatrixXcd Sternheimer::COHSEX()const{
  
  AOBasis dftbasis = _orbitals.SetupDftBasis();
  Vxc_Grid _grid;
  _grid.GridSetup("xxcoarse", _orbitals.QMAtoms(), dftbasis);

  Index nthreads = OPENMP::getMaxThreads();

  //Evaluate static screened Exchange part by integration over r

  Eigen::MatrixXcd Sigma_s =
      Eigen::MatrixXcd::Zero(_density_Matrix.rows(), _density_Matrix.cols());
  for (Index i = 0; i < _grid.getBoxesSize(); ++i) {
    const GridBox& box = _grid[i];
    if (!box.Matrixsize()) {
      continue;
    }

    std::vector<Eigen::MatrixXcd> Sigma_s_thread = std::vector<Eigen::MatrixXcd>(
        nthreads, Eigen::MatrixXcd::Zero(box.Matrixsize(), box.Matrixsize()));

    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

  #pragma omp parallel for schedule(guided)
    for (Index p = 0; p < box.size(); p++) {
      Eigen::VectorXd ao = box.CalcAOValues(points[p]);
      Eigen::MatrixXd S = ao * ao.transpose();
      const double weight = weights[p];

      // Evaluate bare and screend coulomb potential at point to evaluate the
      // correlation screened Coulomb potential (W_c = W-v). This is evaluated
      // in DeltaNsc

      Eigen::MatrixXcd W_s = ScreenedCoulomb(points[p], 0);

      W_s+=-S*W_s;

      Sigma_s_thread[OPENMP::getThreadId()] +=
          weight * S * box.ReadFromBigMatrix(W_s);
          
    }

    Eigen::MatrixXcd Sigma_s_box = std::accumulate(
        Sigma_s_thread.begin(), Sigma_s_thread.end(),
        Eigen::MatrixXcd::Zero(box.Matrixsize(), box.Matrixsize()).eval());
    box.AddtoBigMatrix(Sigma_s, Sigma_s_box);
  }

  //Calculate the static coulomb hole part

  AOBasis auxbasis = _orbitals.SetupAuxBasis();

  AOCoulomb AOC;

  AOC.Fill(dftbasis);

  Eigen::MatrixXcd coulombmatrix = AOC.Matrix();
  Eigen::MatrixXcd DeltaN = DeltaNSC(0, coulombmatrix);
  Eigen::MatrixXcd HartreeInt = _eris.ContractRightIndecesWithMatrix(DeltaN);
  Eigen::MatrixXcd FxcInt =
      Fxc(DeltaN);  // Vxcpot.IntegrateFXC(_density_Matrix, DeltaN);
  // Eigen::MatrixXcd DeltaV = coulombmatrix + HartreeInt + FxcInt;
  Eigen::MatrixXcd DeltaV =
      HartreeInt + FxcInt + coulombmatrix;

  Eigen::MatrixXcd Sigma_c=1/2*DeltaV;

  return Sigma_c+Sigma_s;
}

}  // namespace xtp

}  // namespace votca