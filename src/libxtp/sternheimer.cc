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

#include <fstream>
#include <math.h>
#include <votca/tools/property.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/padeapprox.h>
#include <votca/xtp/multishift.h>
#include <votca/xtp/sternheimer.h>
#include <eigen3/Eigen/src/Core/Matrix.h>

#include "votca/xtp/ERIs.h"

namespace votca {
namespace xtp {

void Sternheimer::Initialize() {

  this->_num_occ_lvls = _orbitals.getNumberOfAlphaElectrons();
  this->_basis_size = _orbitals.getBasisSetSize();
  this->_overlap_Matrix = OverlapMatrix();
  this->_density_Matrix = _orbitals.DensityMatrixGroundState();
  this->_mo_coefficients = _orbitals.MOCoefficients();
  this->_mo_energies = _orbitals.MOEnergies();
  this->_inverse_overlap = _overlap_Matrix.inverse();
  this->_Hamiltonian_Matrix = Hamiltonian();
}    

void Sternheimer::initializeMultishift(int size) { 
  _multishift.setMatrixSize(size);
}
void Sternheimer::initializePade(int size) {
  _pade.initialize(size);
}

Eigen::MatrixXcd Sternheimer::OverlapMatrix() {

  AOBasis basis = _orbitals.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  return overlap.Matrix();
}

Eigen::MatrixXcd Sternheimer::DensityMatrix() {
  return _orbitals.DensityMatrixGroundState();
}

Eigen::MatrixXcd Sternheimer::Hamiltonian() {

  const Eigen::MatrixXcd& mo_coefficients = _orbitals.MOCoefficients();
  const Eigen::MatrixXd& mo_energies = _orbitals.MOEnergies().asDiagonal();
  return _overlap_Matrix * mo_coefficients * mo_energies *
                      mo_coefficients.transpose() * _overlap_Matrix;
}

Eigen::MatrixXcd Sternheimer::CoulombMatrix() {

  AOBasis basis = _orbitals.SetupDftBasis();
  AOCoulomb coulomb;
  coulomb.Fill(basis);
  return coulomb.Matrix();
}

Eigen::MatrixXcd Sternheimer::SternheimerLHS(const Eigen::MatrixXcd& hamiltonian,
                                             const Eigen::MatrixXcd& inverse_overlap,
                                             double eps,
                                             std::complex<double> omega,
                                             bool pm) const{

  std::complex<double> temp=eps+omega;
  if (pm != true) {
    temp = (eps - omega);
  }
  return (inverse_overlap*hamiltonian - temp*Eigen::MatrixXcd::Identity(_basis_size,_basis_size));
}

Eigen::VectorXcd Sternheimer::SternheimerRHS(const Eigen::MatrixXcd& inverse_overlap,
                                             const Eigen::MatrixXcd& density,
                                             const Eigen::MatrixXcd& pertubation,
                                             const Eigen::VectorXcd& coeff) const{
  return -1*(inverse_overlap-density)*pertubation*coeff;
}

std::vector<Eigen::MatrixXcd> Sternheimer::DeltaNOneShot(
    std::vector<std::complex<double>> w,const Eigen::MatrixXd& pertubation) const{
  
  double alpha = 2 * (_mo_energies(_num_occ_lvls) - _mo_energies(2));

  std::vector<Eigen::MatrixXcd> solution_p;
  std::vector<Eigen::MatrixXcd> solution_m;

  Eigen::MatrixXcd H_new;
  Eigen::MatrixXcd LHS_P;
  Eigen::MatrixXcd LHS_M;

  Eigen::VectorXcd RHS;
  
  Multishift::MultiShiftResult result;
  
  for (int v = 0; v < _num_occ_lvls; v++) {

    RHS = SternheimerRHS(_inverse_overlap, _density_Matrix, pertubation, _mo_coefficients.col(v));
    
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
        LHS_M = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w[i], false);
        
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

Eigen::MatrixXcd Sternheimer::DeltaNSelfConsistent(std::complex<double> w,
                                                                const Eigen::MatrixXd& initGuess) const{
    Eigen::MatrixXcd solution_p;
    Eigen::MatrixXcd solution_m;
    
    Eigen::MatrixXd pertubation;
    Eigen::MatrixXcd delta_n;
    AOBasis dftbasis = _orbitals.SetupDftBasis();
    ERIs eris;
    eris.Initialize(dftbasis, dftbasis);
    
    for(int n = 0; n<100; n++){
        for (int v = 0; v < _num_occ_lvls; v++) {

              Eigen::MatrixXcd RHS = SternheimerRHS(_inverse_overlap, _density_Matrix, pertubation, _mo_coefficients.col(v));
              if (w == 0.0) {
                //H_new = _H + alpha * _S * _p.transpose();
                Eigen::MatrixXcd LHS_P = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, true);
                solution_p.col(v)=LHS_P.colPivHouseholderQr().solve(RHS);
              } else {
                Eigen::MatrixXcd LHS_P = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, true);
                solution_p.col(v) = LHS_P.colPivHouseholderQr().solve(RHS);
                }

              if (w == 0.0) {
                //H_new = _H + alpha * _S * _p.transpose();
                Eigen::MatrixXcd LHS_M = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, false);
                solution_m.col(v) = LHS_M.colPivHouseholderQr().solve(RHS);
              }else {
                Eigen::MatrixXcd LHS_M = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, false);
                solution_m.col(v) = LHS_M.colPivHouseholderQr().solve(RHS);
              }
            }
            delta_n =
            2 * _mo_coefficients * solution_p.transpose() +
            2 * _mo_coefficients * solution_m.transpose();

            Eigen::MatrixXd update=eris.ContractRightIndecesWithMatrix(delta_n);
            if(update.norm(),1e-5){
                std::cout<<"Converged at interation n="<<n<<std::endl;
                return delta_n;
            }
                pertubation=initGuess+update;
            }
        std::cout<<"Not Converged"<<std::endl;
        return delta_n;    
    }

std::vector<Eigen::MatrixXcd> Sternheimer::Polarisability(
    const std::vector<std::complex<double>>& grid_w,
    const std::vector<std::complex<double>>& w) {

  initializePade(3);

  AOBasis basis = _orbitals.SetupDftBasis();
  AODipole dipole;
  dipole.Fill(basis);

  std::cout << "Dipole integral complete" << std::endl;
  
  std::vector<Eigen::MatrixXcd> Polar;
  
  for(int n=0; n<grid_w.size(); n++){
      Polar.push_back(Eigen::MatrixXcd::Zero(3, 3));
      for (int i = 0; i < 3; i++) {
        Eigen::MatrixXcd delta_n = DeltaNSelfConsistent(grid_w[n], dipole.Matrix()[i]);
        for (int j = i; j < 3; j++) {
          Polar[n](i, j)+=(delta_n.cwiseProduct(dipole.Matrix()[j])).sum();
        }
      }    
      for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
          Polar[n](j, i)=conj(Polar[n](i, j));
       }
    }
  
      _pade.addPoint(grid_w[n], Polar[n]);
      _pade.addPoint(conj(grid_w[n]), Polar[n].adjoint());
      std::cout<<"Done with w="<<w[n]<<std::endl;
    }
  std::vector<Eigen::MatrixXcd> Polar_pade;
  
  for (std::complex<double> w:w){
    std::cout<<"Calculated Point number"<<w<<std::endl;
    Polar_pade.push_back(_pade.evaluatePoint(w));
  }
  return Polar_pade;
}
}  // namespace xtp
}  // namespace votca

