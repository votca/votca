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

void SternheimerW::Initialize() {

  this->_num_occ_lvls = _orbitals.getNumberOfAlphaElectrons();
  this->_basis_size = _orbitals.getBasisSetSize();
  this->_Hamiltonian_Matrix = Hamiltonian();
  this->_density_Matrix = _orbitals.DensityMatrixGroundState();
  this->_overlap_Matrix = OverlapMatrix();
  this->_mo_coefficients = _orbitals.MOCoefficients();
  this->_mo_energies = _orbitals.MOEnergies();
  this->_inverse_overlap = _overlap_Matrix.inverse();
}

void SternheimerW::initializeMultishift(int size) { 
  _multishift.setMatrixSize(size);
}
void SternheimerW::initializePade(int size) {
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
                                              const Eigen::MatrixXcd& inverse_overlap,
                                              double eps,
                                              std::complex<double> omega,
                                              bool pm) const{

  // distinguish between +w and -w
  std::complex<double> temp=eps+omega;
  if (pm != true) {
    temp = (eps - omega);
  }
  return (inverse_overlap*hamiltonian - temp*Eigen::MatrixXcd::Identity(_basis_size,_basis_size));
}

Eigen::VectorXcd SternheimerW::SternheimerRHS(const Eigen::MatrixXcd& inverse_overlap,
                                              const Eigen::MatrixXcd& density,
                                              const Eigen::MatrixXcd& pertubation,
                                              const Eigen::VectorXcd& coeff) const{

  return -1*(inverse_overlap-density)*pertubation*coeff;
}

Eigen::VectorXcd SternheimerW::SternheimerSolve(const Eigen::MatrixXcd& LHS,
                                                const Eigen::VectorXcd& RHS){
  return LHS.colPivHouseholderQr().solve(RHS);
}

std::vector<Eigen::MatrixXcd> SternheimerW::DeltaNOneShot(
    std::vector<std::complex<double>> w, Eigen::Vector3d r) const{

  Eigen::MatrixXcd V = CoulombMatrix(r);
  
  double alpha = 2 * (_mo_energies(_num_occ_lvls) - _mo_energies(2));

  std::vector<Eigen::MatrixXcd> solution_p;
  std::vector<Eigen::MatrixXcd> solution_m;

  Eigen::MatrixXcd H_new;
  Eigen::MatrixXcd LHS_P;
  Eigen::MatrixXcd LHS_M;

  Eigen::VectorXcd RHS;// = Eigen::Xcd::Zero(_basis_size, _basis_size);

  Multishift::MultiShiftResult result;
  
  for (int v = 0; v < _num_occ_lvls; v++) {

    RHS = SternheimerRHS(_inverse_overlap, _density_Matrix, V, _mo_coefficients.col(v));
    
    for (int i = 0; i < w.size(); i++) {

      if (v == 0) {
        solution_p.push_back(Eigen::MatrixXcd::Zero(_basis_size, _basis_size));
      }

      if (i == 0) {
        //H_new = _H + alpha * _S * _p.transpose();
        LHS_P = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w[i], true);
        // std::cout<<"Sternheimer Setup Complete"<<v<<std::endl;
        //result = _multishift.ComplexBiCG(LHS_P, RHS);
        solution_p[i].col(v)=LHS_P.colPivHouseholderQr().solve(RHS);
        //solution_p[i].col(v) = result._x;
      } else {
        LHS_P = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w[i], true);
        solution_p[i].col(v) = LHS_P.colPivHouseholderQr().solve(RHS);
        //solution_p[i].col(v) = _multishift.DoMultishift(LHS_P,RHS,w[i],result);
//        if(((LHS_P+w[i]*Eigen::MatrixXcd::Identity(_basis_size,_basis_size))*solution_p[i].col(v)-RHS).norm()>1e-13){
//            std::cout<<"res_p="<<(LHS_P+w[i]*Eigen::MatrixXcd::Identity(_basis_size,_basis_size))*solution_p[i].col(v)-RHS<<std::endl;
//        }
      }
    }
    for (int i = 0; i < w.size(); i++) {

      if (v == 0) {
        solution_m.push_back(Eigen::MatrixXcd::Zero(_basis_size, _basis_size));
      }

      if (i == 0) {
        //H_new = _H + alpha * _S * _p.transpose();
        LHS_M = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w[i], false);
        //result = _multishift.ComplexBiCG(LHS_M, RHS);
        solution_m[i].col(v)=LHS_M.colPivHouseholderQr().solve(RHS);
        //solution_m[i].col(v) = result._x;
      }else {
        //LHS_M = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w[i], false);
        
        solution_m[i].col(v) = LHS_M.colPivHouseholderQr().solve(RHS);
        //solution_m[i].col(v) = _multishift.DoMultishift(LHS_M,RHS,-w[i],result);
        //if(((LHS_M+w[i]*Eigen::MatrixXcd::Identity(_basis_size,_basis_size))*solution_m[i].col(v)-RHS).norm()>1e-13){
        //    std::cout<<"res_m="<<(LHS_M-w[i]*Eigen::MatrixXcd::Identity(_basis_size,_basis_size))*solution_m[i].col(v)-RHS<<std::endl;
        //}
      }
    }
  }

  std::vector<Eigen::MatrixXcd> delta_n;
  for (int m = 0; m < w.size(); m++) {

    delta_n.push_back(Eigen::MatrixXcd::Zero(_basis_size, _basis_size));

    delta_n[m] +=
    2 * _mo_coefficients * solution_p[m].transpose() +
    2 * _mo_coefficients * solution_m[m].transpose();
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

  std::cout << std::endl<< "gridsize= " << numint.getGridSize() << std::endl;
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
    V = CoulombMatrix(*gridpoints[r].second);

    for (int v = 0; v < num_occ_lvls; v++) {

      RHS = SternheimerRHS(S, p, V, mo_coefficients.col(v));

      delta_c_p.col(v) = SternheimerSolve(LHS_P[v], RHS);
      delta_c_m.col(v) = SternheimerSolve(LHS_M[v], RHS);
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
    const std::vector<std::complex<double>>& grid_w,
    const std::vector<std::complex<double>>& w, std::string gridtype) {

  //initializePade(3);
  initializeMultishift(_orbitals.getBasisSetSize());

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  NumericalIntegration numint;
  numint.GridSetup(gridtype, _orbitals.QMAtoms(),
                   dftbasis);  // For now use medium grid
  std::vector<std::pair<double, const Eigen::Vector3d*>> grid =
      numint.getGridpoints();

  std::cout << "gridsize= " << numint.getGridSize() << std::endl;

  AOBasis basis = _orbitals.SetupDftBasis();
  AODipole dipole;
  dipole.Fill(basis);

  std::cout << "Dipole integral complete" << std::endl;
  
  std::vector<Eigen::MatrixXcd> Polar;
  for(const std::complex<double> w: grid_w){
      Polar.push_back(Eigen::MatrixXcd::Zero(3, 3));
  }
  
  int index=0;
  for (const std::pair<double, const Eigen::Vector3d*> point: grid) {  
    std::vector<Eigen::MatrixXcd> delta_n = DeltaNOneShot(grid_w, *point.second);
    Eigen::Vector3d gridcont = *point.second;
    for (int o = 0; o < grid_w.size(); o++) {
      for (int i = 0; i < 3; i++) {
        double r_x = gridcont(i);
        for (int j = i; j < 3; j++) {
            Polar[o](i, j)+=point.first*r_x*(delta_n[o].cwiseProduct(dipole.Matrix()[j])).sum();
        }
      }
    }
    if(index%1000==0){
        std::cout<<"Iteration="<<index<<std::endl;
    }
    index++;
    delta_n.clear();
  }
  for (int o = 0; o < grid_w.size(); o++) {
    for (int i = 0; i < 3; i++) {
      for (int j = i+1; j < 3; j++) {
        Polar[o](j, i)=conj(Polar[o](i, j));
      }
    }
  }
  
  for (int i = 0; i < Polar.size(); i++) {
    //_pade.addPoint(grid_w[i], Polar[i].real());
    //_pade.addPoint(conj(grid_w[i]), Polar[i].real());
  }

  //std::vector<Eigen::MatrixXcd> Polar_pade;
  for (std::complex<double> w:w) {
    std::cout<<"Calculated Point number"<<w<<std::endl;
    //Polar_pade.push_back(_pade.evaluatePoint(w));
  }
  return Polar;
}
  std::vector<Eigen::MatrixXcd> SternheimerW::DeltaNSelfConsistent(std::vector<std::complex<double>>& frequency, Eigen::Vector3d& r) const{
      
      
      std::vector<Eigen::MatrixXcd> delta_V;
      
      std::vector<Eigen::MatrixXcd> solution_p;
      std::vector<Eigen::MatrixXcd> solution_m;
      std::vector<Eigen::MatrixXcd> delta_n;
      
      int max_iter=1e4;
      double tol=1e-9;
      for(int n=1;n<max_iter;n++){
          for (int i = 0; i< frequency.size(); i++) {
            delta_V.push_back(CoulombMatrix(r));
            solution_p.push_back(Eigen::MatrixXcd::Zero(_basis_size, _basis_size));
            solution_m.push_back(Eigen::MatrixXcd::Zero(_basis_size, _basis_size));
            for (int v = 0; v < _num_occ_lvls; v++){
                Eigen::MatrixXcd RHS = SternheimerRHS(_inverse_overlap, _density_Matrix, delta_V[i], _mo_coefficients.col(v));
                Eigen::MatrixXcd LHS_P = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), frequency[i], true);
                Eigen::MatrixXcd LHS_M = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), frequency[i], false);
                solution_p[i].col(v) = LHS_P.colPivHouseholderQr().solve(RHS);
                solution_m[i].col(v) = LHS_M.colPivHouseholderQr().solve(RHS);
            }
            delta_n.push_back(Eigen::MatrixXcd::Zero(_basis_size, _basis_size));

            delta_n[i] +=
            2 * _mo_coefficients * solution_p[i].transpose() +
            2 * _mo_coefficients * solution_m[i].transpose();
            delta_V[i] += delta_n[i]*_overlap_Matrix*CoulombMatrix(r);
            if((delta_n[i]*_overlap_Matrix*CoulombMatrix(r)).norm()<tol){
                
                return delta_n;
                
            }
          }
      }         
  }
}  // namespace xtp
}  // namespace votca