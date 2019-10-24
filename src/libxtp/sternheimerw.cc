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

#include <ctime>
#include <fstream>
#include <math.h>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/multishift.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/padeapprox.h>
#include <votca/xtp/sternheimerw.h>

#include "votca/xtp/numerical_integrations.h"

namespace votca {
namespace xtp {

void SternheimerW::Setup() {

  this->_num_occ_lvls = _orbitals.getNumberOfAlphaElectrons();
  this->_basis_size = _orbitals.getBasisSetSize();
  this->_H = Hamiltonian();
  this->_p = _orbitals.DensityMatrixGroundState();
  this->_S = OverlapMatrix();
  this->_mo_coefficients = _orbitals.MOCoefficients();
  this->_mo_energies = _orbitals.MOEnergies();
}

void SternheimerW::initializeMultishift(int basis_size) {

  Multishift multi;
  this->_multishift = multi;

  _multishift.setBasisSize(basis_size);
}
void SternheimerW::initializePade(int size) {

  PadeApprox pade;
  this->_pade = pade;

  _pade.initialize(size);
}

Eigen::MatrixXd SternheimerW::OverlapMatrix() const{

  AOBasis basis = _orbitals.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  return overlap.Matrix();
}

Eigen::MatrixXd SternheimerW::DensityMatrix() const{

  return _orbitals.DensityMatrixGroundState();
}

Eigen::MatrixXd SternheimerW::Hamiltonian() const{

  const Eigen::MatrixXd& mo_coefficients = _orbitals.MOCoefficients();
  const Eigen::MatrixXd& mo_energies = _orbitals.MOEnergies().asDiagonal();
  const Eigen::MatrixXd overlap = OverlapMatrix();
  Eigen::MatrixXd H = overlap * mo_coefficients * mo_energies *
                      mo_coefficients.transpose() * overlap;

  return H;
}
Eigen::MatrixXcd SternheimerW::CoulombMatrix(Eigen::Vector3d gridpoint) const{

  AOBasis basis = _orbitals.SetupDftBasis();
  AOESP aoesp;
  aoesp.setPosition(gridpoint);
  aoesp.Fill(basis);
  return aoesp.Matrix();
}

Eigen::MatrixXcd SternheimerW::SternheimerLHS(const Eigen::MatrixXcd& hamiltonian,
                                              const Eigen::MatrixXcd& overlap,
                                              double eps,
                                              std::complex<double> omega,
                                              bool pm) const{

  Eigen::MatrixXcd S;
  // distinguish between +w and -w
  if (pm == true) {
    S = (eps + omega) * overlap;
  } else {
    S = (eps - omega) * overlap;
  }
  return (hamiltonian - S);
}

Eigen::VectorXcd SternheimerW::SternheimerRHS(const Eigen::MatrixXcd& overlap,
                                              const Eigen::MatrixXcd& density,
                                              const Eigen::MatrixXcd& pertubation,
                                              const Eigen::VectorXcd& coeff) const{

  Eigen::MatrixXcd I =
      Eigen::MatrixXcd::Identity(overlap.rows(), overlap.cols());
  // Calc RHS
  Eigen::MatrixXcd t = -1 * (I - overlap * density);

  t = t * pertubation;

  Eigen::VectorXcd M = t * coeff;

  return M;
}

Eigen::VectorXcd SternheimerW::SternheimerSolve(const Eigen::MatrixXcd& LHS,
                                                const Eigen::VectorXcd& RHS){

  // Eigen::VectorXcd x=_multishift.CBiCG(LHS, RHS);

  return LHS.colPivHouseholderQr().solve(RHS);;
}

std::vector<Eigen::MatrixXcd> SternheimerW::DeltaNOSOP(
    std::vector<std::complex<double>> w, Eigen::Vector3d r) const{

  Eigen::MatrixXcd V = CoulombMatrix(r);

  double alpha = 2 * (_mo_energies(_num_occ_lvls) - _mo_energies(2));

  std::vector<Eigen::MatrixXcd> delta_c_p;
  std::vector<Eigen::MatrixXcd> delta_c_m;

  Eigen::MatrixXcd H_new;
  Eigen::MatrixXcd LHS_P;
  Eigen::MatrixXcd LHS_M;

  Eigen::VectorXcd RHS;// = Eigen::Xcd::Zero(_basis_size, _basis_size);

  std::vector<Eigen::MatrixXcd> delta_n;

  Multishift::MultiShiftResult result;
  
  for (int v = 0; v < _num_occ_lvls; v++) {

    RHS = SternheimerRHS(_S, _p, V, _mo_coefficients.col(v));

    for (int i = 0; i < w.size(); i++) {

      if (v == 0) {
        delta_c_p.push_back(Eigen::MatrixXcd::Zero(_basis_size, _basis_size));
      }

      if (i == 0) {
        H_new = _H + alpha * _S * _p.transpose();
        LHS_P = SternheimerLHS(H_new, _S, _mo_energies(v), w.at(i), true);
        // std::cout<<"Sternheimer Setup Complete"<<v<<std::endl;
        result = _multishift.CBiCG(LHS_P, RHS);
        delta_c_p.at(i).col(v) = result._x;
      } else {
        LHS_P = SternheimerLHS(_H, _S, _mo_energies(v), w.at(i), true);

        delta_c_p.at(i).col(v) = _multishift.DoMultishift(LHS_P, RHS, w.at(i), result);
      }
    }
    for (int i = 0; i < w.size(); i++) {

      if (v == 0) {
        delta_c_m.push_back(Eigen::MatrixXcd::Zero(_basis_size, _basis_size));
      }

      if (i == 0) {
        H_new = _H + alpha * _S * _p.transpose();
        LHS_M = SternheimerLHS(H_new, _S, _mo_energies(v), w.at(i), false);
        result = _multishift.CBiCG(LHS_M, RHS);
        delta_c_m.at(i).col(v) = result._x;
      }else {
        LHS_M = SternheimerLHS(_H, _S, _mo_energies(v), w.at(i), false);

        delta_c_m.at(i).col(v) = _multishift.DoMultishift(LHS_M, RHS, w.at(i), result);
      }
    }
  }

  for (int m = 0; m < w.size(); m++) {

    delta_n.push_back(Eigen::MatrixXcd::Zero(_basis_size, _basis_size));

    for (int v = 0; v < _basis_size; v++) {

      for (int i = 0; i < _basis_size; i++) {

        for (int j = 0; j < _basis_size; j++) {

          delta_n.at(m)(i, j) =
              delta_n.at(m)(i, j) +
              2 * _mo_coefficients(i, v) * delta_c_p.at(m)(j, v) +
              2 * _mo_coefficients(i, v) * delta_c_m.at(m)(j, v);
        }
      }
    }
  }
  return delta_n;
}

std::vector<Eigen::MatrixXcd> SternheimerW::DeltaNOS(std::complex<double> w,
                                                     std::string gridtype) {

  // Setting up grid for evaluation
  AOBasis dftbasis = _orbitals.SetupDftBasis();
  NumericalIntegration numint;
  numint.GridSetup(gridtype, _orbitals.QMAtoms(),
                   dftbasis);  // For now use medium grid

  std::vector<std::pair<double, const Eigen::Vector3d*>> gridpoints =
      numint.getGridpoints();

  std::cout << "gridsize= " << numint.getGridSize() << std::endl;
  // Setting up constants
  const int& num_occ_lvls = _orbitals.getNumberOfAlphaElectrons();
  const int& basis_size = _orbitals.getBasisSetSize();

  // Setting up constant matrices and vectors needed for the Sternheimer
  // equation
  const Eigen::MatrixXd H = Hamiltonian();
  const Eigen::MatrixXd S = OverlapMatrix();
  const Eigen::MatrixXd p = DensityMatrix();
  Eigen::MatrixXcd V;

  const Eigen::MatrixXd& mo_coefficients = _orbitals.MOCoefficients();
  const Eigen::VectorXd& mo_energies = _orbitals.MOEnergies();

  // Setting up container for solutions in each gridpoint

  Eigen::MatrixXcd delta_c_p = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::MatrixXcd delta_c_m = Eigen::MatrixXcd::Zero(basis_size, basis_size);

  // Setting up container for LHS and RHS of the system
  std::vector<Eigen::MatrixXcd> LHS_P;
  std::vector<Eigen::MatrixXcd> LHS_M;
  Eigen::VectorXcd RHS;

  // H_new is used if omega=0 to shift the eigenvalues away from 0
  Eigen::MatrixXcd H_new;
  double alpha = 2 * (mo_energies(num_occ_lvls) - mo_energies(2));

  //double pi = tools::conv::Pi;
  // Initilazing solution vector and container for solution

  std::vector<Eigen::MatrixXcd> result;
  Eigen::MatrixXcd delta_n = Eigen::MatrixXcd::Zero(basis_size, basis_size);

  // Iterating over all occupied states and saving the LHS
  for (int v = 0; v < num_occ_lvls; v++) {

    // Setting up LHS since it is not dependent on r
    if (w == 0.0) {
      H_new = H + alpha * S * p.transpose();
      LHS_P.push_back(SternheimerLHS(H_new, S, mo_energies(v), w, true));
      LHS_M.push_back(SternheimerLHS(H_new, S, mo_energies(v), w, false));
    } else {
      LHS_P.push_back(SternheimerLHS(H, S, mo_energies(v), w, true));
      LHS_M.push_back(SternheimerLHS(H, S, mo_energies(v), w, false));
    }
  }
  // Iterate over all spacial gridpoints
  for (int r = 0; r < gridpoints.size(); r++) {

    // std::cout<<(*gridpoints[r])<<std::endl;
    V = CoulombMatrix(*gridpoints.at(r).second);

    for (int v = 0; v < num_occ_lvls; v++) {

      RHS = SternheimerRHS(S, p, V, mo_coefficients.col(v));

      delta_c_p.col(v) = SternheimerSolve(LHS_P.at(v), RHS);
      delta_c_m.col(v) = SternheimerSolve(LHS_M.at(v), RHS);
    }

    for (int v = 0; v < basis_size; v++) {
      for (int i = 0; i < basis_size; i++) {
        for (int j = 0; j < basis_size; j++) {
          delta_n(i, j) = delta_n(i, j) +
                          2 * mo_coefficients(i, v) * delta_c_p(j, v) +
                          2 * mo_coefficients(i, v) * delta_c_m(j, v);
        }
      }
    }

    result.push_back(delta_n);
    delta_n = Eigen::MatrixXcd::Zero(basis_size, basis_size);

    if (r % 100 == 0) {

      std::cout << "I'm in iteration " << r << std::endl;
    }
  }

  return result;
}

std::vector<Eigen::MatrixXcd> SternheimerW::Polarisability(
    std::vector<std::complex<double>> grid_w,
    std::vector<std::complex<double>> w, std::string gridtype) {

  initializePade(3);
  initializeMultishift(_orbitals.getBasisSetSize());

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  NumericalIntegration numint;
  numint.GridSetup(gridtype, _orbitals.QMAtoms(),
                   dftbasis);  // For now use medium grid
  const int& basis_size = _orbitals.getBasisSetSize();
  std::vector<std::pair<double, const Eigen::Vector3d*>> grid =
      numint.getGridpoints();

  std::cout << "gridsize= " << numint.getGridSize() << std::endl;

  std::vector<Eigen::MatrixXcd> Polar;
  std::vector<Eigen::MatrixXcd> Polar_pade;

  std::vector<Eigen::MatrixXcd> delta_n;

  double r_x;
  Eigen::Vector3d gridcont;

  AOBasis basis = _orbitals.SetupDftBasis();
  AODipole dipole;
  dipole.Fill(basis);

  std::cout << "Dipole integral complete" << std::endl;

  for (int r = 0; r < grid.size(); r++) {

    delta_n = DeltaNOSOP(grid_w, *grid.at(r).second);

    _pade.clear();

    for (int o = 0; o < grid_w.size(); o++) {
      if (r == 0) Polar.push_back(Eigen::MatrixXcd::Zero(3, 3));
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          for (int n = 0; n < basis_size; n++) {
            for (int m = 0; m < basis_size; m++) {
              gridcont = *grid.at(r).second;
              r_x = gridcont(i);
              Polar.at(o)(i, j) = Polar.at(o)(i, j) +
                                  delta_n.at(o)(n, m) * r_x * grid.at(r).first *
                                      dipole.Matrix()[j](n, m);
            }
          }
        }
      }
    }
    if (r % 100 == 0) {
      std::cout << "I'm in iteration " << r << std::endl;
    }
    delta_n.clear();
  }
  for (int i = 0; i < Polar.size(); i++) {
    // std::cout<<"Complex Polar at
    // w="<<grid_w.at(i)<<std::endl<<Polar.at(i)<<std::endl;
    _pade.addPoint(grid_w.at(i), Polar.at(i));
    // std::cout<<"Complex Polar* at
    // w*="<<conj(grid_w.at(i))<<std::endl<<Polar.at(i).adjoint()<<std::endl;
    _pade.addPoint(conj(grid_w.at(i)), Polar.at(i).adjoint());
  }

  for (int i = 0; i < w.size(); i++) {
    Polar_pade.push_back(_pade.evaluate(w.at(i)));
  }
  return Polar_pade;
}

void SternheimerW::testPade() {

  PadeApprox pade;

  pade.test();
}

bool SternheimerW::evaluate() { return true; }
}  // namespace xtp
}  // namespace votca