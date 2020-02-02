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
#include <fstream>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aomatrix3d.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/multishift.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/padeapprox.h>
#include <votca/xtp/sternheimer.h>
#include <votca/xtp/vxc_potential.h>
#include <votca/xtp/vxc_grid.h>

namespace votca {
namespace xtp {

void Sternheimer::Initialize() {

  this->_num_occ_lvls = _orbitals.getNumberOfAlphaElectrons();
  this->_basis_size = _orbitals.getBasisSetSize();
  this->_overlap_Matrix = OverlapMatrix();
  this->_density_Matrix = _orbitals.DensityMatrixGroundState();
  this->_mo_coefficients = _orbitals.MOs().eigenvectors();
  this->_mo_energies = _orbitals.MOs().eigenvalues();
  this->_inverse_overlap = _overlap_Matrix.inverse();
  this->_Hamiltonian_Matrix = Hamiltonian();
}

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

Eigen::MatrixXcd Sternheimer::DeltaNOneShot(
    std::complex<double> w, const Eigen::MatrixXcd& perturbation) const {
  Eigen::MatrixXcd solution_p =
      Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);
  Eigen::MatrixXcd solution_m =
      Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);

  double e_field = 1E-1;
  double mixing = 0.5;
  Index maxerrorindex = 5;

  std::vector<Eigen::MatrixXcd> perturbationVectorInput;
  std::vector<Eigen::MatrixXcd> perturbationVectoroutput;
  perturbationVectorInput.push_back(-e_field * perturbation);
  Eigen::MatrixXcd perturbationUsed = (-e_field * perturbation);
  // Eigen::MatrixXcd delta_n_in_new = Eigen::MatrixXcd::Zero(_basis_size,
  // _basis_size); Eigen::MatrixXcd delta_n_in_old =
  // Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  Eigen::MatrixXcd delta_n_out_new =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  Eigen::MatrixXcd delta_n_out_old =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  Eigen::MatrixXcd delta_n_step_one =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  Eigen::MatrixXcd V_ext = -e_field * perturbation;

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  AOBasis auxbasis = _orbitals.SetupAuxBasis();
  ERIs eris;
  eris.Initialize(dftbasis, auxbasis);
  // eris.Initialize(dftbasis);

  // ADIIS adiis;
  // DIIS diis;
  // diis.setHistLength(10);
  // diis.Update(0,10e9*(-e_field * perturbation));

Vxc_Grid grid;
grid.GridSetup("fine",_orbitals.QMAtoms(),dftbasis);
Vxc_Potential<Vxc_Grid>Vxcpot(grid);
std::cout<<_orbitals.getXCFunctionalName()<<std::endl;
Vxcpot.setXCfunctional(_orbitals.getXCFunctionalName());

  double diff = 10000;

  for (Index n = 0; n < 200; n++) {

    for (Index v = 0; v < _num_occ_lvls; v++) {
      // std::cout<<"Loop v="<<v<<std::endl;
      Eigen::MatrixXcd RHS =
          SternheimerRHS(_inverse_overlap, _density_Matrix, perturbationUsed,
                         _mo_coefficients.col(v));

      Eigen::MatrixXcd LHS_P = SternheimerLHS(
          _Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, true);
      solution_p.col(v) = LHS_P.colPivHouseholderQr().solve(RHS);
      Eigen::MatrixXcd LHS_M = SternheimerLHS(
          _Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, false);
      solution_m.col(v) = LHS_M.colPivHouseholderQr().solve(RHS);
    }
    // std::cout<<"Solving Sternheimer complete"<<std::endl;
    delta_n_out_old = delta_n_out_new;
    delta_n_out_new =
        2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
            solution_p.transpose() +
        2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
            solution_m.transpose();

    Eigen::MatrixXcd delta = (delta_n_out_new - delta_n_out_old);
    diff = delta.squaredNorm();

    // if(diff<10e-9){
    //   std::cout<<"Converged after " << n+1 <<" iteration."<<std::endl;
    //   //throw std::exception();
    //   return delta_n_out_new;
    // }

    Eigen::MatrixXcd contract =
        eris.ContractRightIndecesWithMatrix(delta_n_out_new);

    Eigen::MatrixXcd FxcInt = Vxcpot.IntegrateFXC(_density_Matrix,delta_n_out_new);

    std::cout<<"Norm of Fxc"<<FxcInt.norm()<<std::endl;

    if (perturbationVectoroutput.size() > 4) {
      perturbationVectoroutput.erase(perturbationVectoroutput.begin());
    }

    perturbationVectoroutput.push_back((-e_field * perturbation) + contract +FxcInt);

    diff = (perturbationVectorInput.back() - perturbationVectoroutput.back())
               .squaredNorm();
    std::cout << n << " " << diff << std::endl;
    if (diff < 10e-9) {
      std::cout << "Converged after " << n + 1 << " iteration. TEST" << std::endl;
      // throw std::exception();
      return delta_n_out_new;
    }

    // perturbationUsed =
    // mixing*(perturbationVector.at(perturbationVector.size()-2))+(1-mixing)*perturbationVector.at(perturbationVector.size()-1);
    if (n == 0) {
      perturbationUsed =
          mixing * perturbationVectoroutput.back() +
          (1 - mixing) * perturbationVectorInput
                             .back();  // at(perturbationVectorInput.size()-1);
      perturbationVectorInput.push_back(perturbationUsed);
    } else {

      // perturbationUsed = (BroydenMixing(perturbationVectorInput,
      // perturbationVectoroutput, 0.5));
      perturbationUsed = (NPAndersonMixing(perturbationVectorInput,
                                           perturbationVectoroutput, 0.5));
      // std::cout<<"test"<<std::endl;
      if (perturbationVectorInput.size() > 4) {
        perturbationVectorInput.erase(perturbationVectorInput.begin());
      }
      perturbationVectorInput.push_back(perturbationUsed);
      // std::cout<<"test2"<<std::endl;
    }
  }
  std::cout << "NOT converged diff = " << diff << " The frequency is w = " << w
            << std::endl;
  // std::cout<<"Returning delta n from first iteration"<<std::endl;
  return delta_n_step_one;
}

Eigen::MatrixXcd Sternheimer::AndersonMixing(Eigen::MatrixXcd inNew,
                                             Eigen::MatrixXcd inOld,
                                             Eigen::MatrixXcd outNew,
                                             Eigen::MatrixXcd outOld,
                                             double alpha) const {

  std::complex<double> beta =
      (outNew - inNew).cwiseProduct(outNew - inNew - outOld + inOld).sum() /
      ((outNew - inNew).cwiseProduct((outOld - inOld))).sum();

  Eigen::MatrixXcd nIn = beta * inOld + (1 - beta) * inNew;
  Eigen::MatrixXcd nOut = beta * outOld + (1 - beta) * outNew;

  // std::cout<<"beta = "<<beta<<std::endl;

  return alpha * nOut + (1 - alpha) * nIn;
}

Eigen::MatrixXcd Sternheimer::NPAndersonMixing(
    std::vector<Eigen::MatrixXcd>& Input, std::vector<Eigen::MatrixXcd>& Output,
    double alpha) const {

  // std::cout<<"Started Anderson Mixing with size "<<Input.size()<<" Output
  // size = "<<Output.size()<<std::endl;

  // std::vector<Eigen::VectorXcd> Input_vec;
  // std::vector<Eigen::VectorXcd> Output_vec;

  // for(Index i=0; i<Input.size(); i++){

  //   Eigen::Map<Eigen::VectorXcd> v1(Input.at(i).data(), Input.at(i).size());
  //   Eigen::Map<Eigen::VectorXcd> v2(Output.at(i).data(),
  //   Output.at(i).size()); Input_vec.push_back(v1); Output_vec.push_back(v2);

  // }

  // Index size = Input_vec.size();

  // Calculating Delta N and saving it for speedup

  Eigen::MatrixXcd DeltaN = Output.back() - Input.back();
  // Eigen::VectorXcd DeltaN= Output_vec.back()-Input.back();

  // std::cout<<"1"<<std::endl;

  // Building Linear System for Coefficients
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(Input.size() - 1, Input.size() - 1);
  Eigen::VectorXd c = Eigen::VectorXd::Zero(Input.size() - 1);

  // std::cout<<"2"<<std::endl;

  for (Index m = 1; m < Input.size(); m++) {

    c(m - 1) = (DeltaN - Output.at(Output.size() - 1 - m) +
                Input.at(Input.size() - 1 - m))
                   .cwiseProduct(DeltaN)
                   .sum()
                   .real();
    // c(m-1)=(DeltaN-Output_vec.at(size-1-m)+Input_vec.at(size-1-m)).dot(DeltaN).real();
    for (Index j = 1; j < Input.size(); j++) {

      A(m - 1, j - 1) =
          (DeltaN - Output.at(Output.size() - 1 - m) +
           Input.at(Input.size() - 1 - m))
              .cwiseProduct(DeltaN - Output.at(Output.size() - 1 - j) +
                            Input.at(Input.size() - 1 - j))
              .sum()
              .real();
      // A(m-1,j-1) =
      // (DeltaN-Output_vec.at(size-1-m)+Input_vec.at(size-1-m)).dot(DeltaN-Output_vec.at(size-1-j)+Input_vec.at(size-1-j));
    }
  }
  // Solving the System to obtain coefficients
  Eigen::VectorXcd coefficients = A.fullPivHouseholderQr().solve(c);

  // Mixing the Potentials

  Eigen::MatrixXcd OutMixed = Output.back();
  Eigen::MatrixXcd InMixed = Input.back();
  // Eigen::VectorXcd OutMixed = Output_vec.back();
  // Eigen::VectorXcd InMixed = Input_vec.back();

  for (Index n = 1; n < Input.size(); n++) {

    OutMixed += coefficients(n - 1) * (Output.at(Output.size() - 1 - n) -
                                       Output.at(Output.size() - 1));
    InMixed += coefficients(n - 1) *
               (Input.at(Input.size() - 1 - n) - Input.at(Input.size() - 1));
    // OutMixed+=coefficients(n-1)*(Output_vec.at(size-1-n)-Output_vec.back();
    // InMixed+=coefficients(n-1)*(Input_vec.at(size-1-n)-Input_vec.back();
  }

  // Returning the linear Mix of Input and Output
  return alpha * OutMixed + (1 - alpha) * InMixed;
}

Eigen::MatrixXcd Sternheimer::BroydenMixing(
    std::vector<Eigen::MatrixXcd> Input, std::vector<Eigen::MatrixXcd> Output,
    double jacobianScaling) const {

  // Approximation of the first inverse Jacobian
  Eigen::MatrixXcd FirstInverseJacobian =
      jacobianScaling * Eigen::MatrixXcd::Identity(_basis_size, _basis_size);

  Index histSize = Input.size();

  Eigen::MatrixXd gamma = Eigen::MatrixXd::Zero(histSize, histSize);
  Eigen::MatrixXd alpha = Eigen::MatrixXd::Zero(histSize - 1, histSize - 1);
  Eigen::MatrixXd beta = Eigen::MatrixXd::Zero(histSize - 1, histSize - 1);
  Eigen::VectorXd weights = Eigen::VectorXd::Zero(histSize);
  // std::cout<<"1"<<std::endl;
  for (Index m = 0; m < histSize; m++) {
    weights(m) = 1 / sqrt(abs((Output.at(m) - Input.at(m))
                                  .cwiseProduct((Output.at(m) - Input.at(m)))
                                  .sum()
                                  .real()));
    // std::cout<<"check
    // "<<((Output.at(m)-Input.at(m)).cwiseProduct(Output.at(m)-Input.at(m)).real().sum())<<std::endl;
  }
  // std::cout<<"2"<<std::endl;
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

      // if(alpha(m,l)!=alpha(m,l)){
      //   std::cout<<"denominator1 = "<<
      //   (Output.at(l+1)-Input.at(l+1)-Output.at(l)+Input.at(l)).norm()
      //   <<std::endl; std::cout<<"denominator2 = "<<
      //   (Output.at(m+1)-Input.at(m+1)-Output.at(m)+Input.at(m)).norm()
      //   <<std::endl; std::cout<<"weights="<<weights(m)<<"
      //   "<<weights(l)<<std::endl;
      // }
    }
  }
  // std::cout<<"3"<<std::endl;
  // if(abs((weights(0)*weights(0)*Eigen::MatrixXd::Identity(histSize-1,histSize-1)+alpha).determinant())<10e-6){
  // std::cout<<"w = "<<weights(0)<<std::endl;
  // std::cout<<"alpha "<<alpha<<std::endl;

  // }
  beta = (weights(0) * weights(0) *
              Eigen::MatrixXd::Identity(histSize - 1, histSize - 1) +
          alpha)
             .inverse();
  // std::cout<<"3.5"<<std::endl;
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
  // std::cout<<"4"<<std::endl;
  Eigen::MatrixXcd BroydenMix =
      Input.back() + FirstInverseJacobian * (Output.back() - Input.back());
  // std::cout<<"5"<<std::endl;
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

std::vector<Eigen::Matrix3cd> Sternheimer::Polarisability(
    double omega_start, double omega_end, Index steps, double imaginary_shift,
    double lorentzian_broadening, Index resolution_output) const {

  std::vector<std::complex<double>> grid_w =
      BuildGrid(omega_start, omega_end, steps, imaginary_shift);

  for (Index i = 0; i < grid_w.size(); i++) {
    std::cout << grid_w[i] << std::endl;
  }

  std::vector<std::complex<double>> frequency_grid = BuildGrid(
      omega_start, omega_end, resolution_output, lorentzian_broadening);

  std::vector<Eigen::Matrix3cd> Polar;
  std::vector<Eigen::Matrix3cd> Polar_pade;

  for (Index i = 0; i < grid_w.size(); i++) {
    Polar.push_back(Eigen::Matrix3cd::Zero());
  }

  PadeApprox pade_1;
  PadeApprox pade_2;
  PadeApprox pade_3;
  PadeApprox pade_4;
  PadeApprox pade_5;
  PadeApprox pade_6;
  pade_1.initialize(4 * grid_w.size());
  pade_2.initialize(4 * grid_w.size());
  pade_3.initialize(4 * grid_w.size());
  pade_4.initialize(4 * grid_w.size());
  pade_5.initialize(4 * grid_w.size());
  pade_6.initialize(4 * grid_w.size());

  AOBasis basis = _orbitals.SetupDftBasis();
  AODipole dipole;
  dipole.Fill(basis);

#pragma omp parallel for
  for (Index n = 0; n < grid_w.size(); n++) {
    for (Index i = 0; i < 3; i++) {
      Eigen::MatrixXcd delta_n = DeltaNOneShot(grid_w[n], dipole.Matrix()[i]);
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
    pade_1.addPoint(grid_w[n], Polar[n](0, 0));
    pade_1.addPoint(conj(grid_w[n]), conj(Polar[n](0, 0)));
    pade_1.addPoint(-grid_w[n], conj(Polar[n](0, 0)));
    pade_1.addPoint(-conj(grid_w[n]), Polar[n](0, 0));

    // pade_2.addPoint(grid_w[n], Polar[n](0,1));
    // pade_2.addPoint(conj(grid_w[n]), conj(Polar[n](0,1)));
    // pade_2.addPoint(-grid_w[n], conj(Polar[n](0, 1)));
    // pade_2.addPoint(-conj(grid_w[n]), Polar[n](0, 1));

    // pade_3.addPoint(grid_w[n], Polar[n](0,2));
    // pade_3.addPoint(conj(grid_w[n]), conj(Polar[n](0,2)));
    // pade_3.addPoint(-grid_w[n], conj(Polar[n](0, 2)));
    // pade_3.addPoint(-conj(grid_w[n]), Polar[n](0, 2));

    pade_4.addPoint(grid_w[n], Polar[n](1, 1));
    pade_4.addPoint(conj(grid_w[n]), conj(Polar[n](1, 1)));
    pade_4.addPoint(-grid_w[n], conj(Polar[n](1, 1)));
    pade_4.addPoint(-conj(grid_w[n]), Polar[n](1, 1));

    // pade_5.addPoint(grid_w[n], Polar[n](1,2));
    // pade_5.addPoint(conj(grid_w[n]), conj(Polar[n](1,2)));
    // pade_5.addPoint(-grid_w[n], conj(Polar[n](1, 2)));
    // pade_5.addPoint(-conj(grid_w[n]), Polar[n](1, 2));

    pade_6.addPoint(grid_w[n], Polar[n](2, 2));
    pade_6.addPoint(conj(grid_w[n]), conj(Polar[n](2, 2)));
    pade_6.addPoint(-grid_w[n], conj(Polar[n](2, 2)));
    pade_6.addPoint(-conj(grid_w[n]), Polar[n](2, 2));
  }

  for (std::complex<double> w : frequency_grid) {
    Polar_pade.push_back(Eigen::Matrix3cd::Zero());
    Polar_pade[Polar_pade.size() - 1](0, 0) = pade_1.evaluatePoint(w);
    // Polar_pade[Polar_pade.size() - 1](0,1)=pade_2.evaluatePoint(w);
    // Polar_pade[Polar_pade.size() - 1](0,2)=pade_3.evaluatePoint(w);
    // Polar_pade[Polar_pade.size() -
    // 1](1,0)=Polar_pade[Polar_pade.size()-1](0,1);
    Polar_pade[Polar_pade.size() - 1](1, 1) = pade_4.evaluatePoint(w);
    // Polar_pade[Polar_pade.size() - 1](1,2)=pade_5.evaluatePoint(w);
    // Polar_pade[Polar_pade.size() -
    // 1](2,0)=Polar_pade[Polar_pade.size()-1](0,2);
    // Polar_pade[Polar_pade.size()
    // - 1](2,1)=Polar_pade[Polar_pade.size()-1](1,2);
    Polar_pade[Polar_pade.size() - 1](2, 2) = pade_6.evaluatePoint(w);
  }
  printIsotropicAverage(Polar_pade, frequency_grid);
  return Polar_pade;
}
void Sternheimer::printIsotropicAverage(
    std::vector<Eigen::Matrix3cd>& polar,
    std::vector<std::complex<double>>& grid) const {
  for (Index i = 0; i < polar.size(); i++) {
    std::cout << real(grid.at(i)) * votca::tools::conv::hrt2ev << " "
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

}  // namespace xtp
}  // namespace votca