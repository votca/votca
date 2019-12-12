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
#include <votca/xtp/aopotential.h>
#include <votca/xtp/aomatrix3d.h>
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
  this->_mo_coefficients = _orbitals.MOs().eigenvectors();
  this->_mo_energies = _orbitals.MOs().eigenvalues();
  this->_inverse_overlap = _overlap_Matrix.inverse();
}

std::vector<std::complex<double>> SternheimerW::BuildGrid(
    double omega_start, double omega_end, int steps,
    double imaginary_shift) const {

  std::vector<std::complex<double>> grid;

  double stepsize = (omega_end - omega_start) / steps;
  const double ev2hrt = votca::tools::conv::ev2hrt;
  std::complex<double> d(1, 0);
  std::complex<double> i(0, 1);

  for (int n = 0; n <= steps; n++) {
    // Converts from input eV to Hartree
    grid.push_back(omega_start * ev2hrt + n * stepsize * ev2hrt +
                   imaginary_shift * i * ev2hrt);
  }

  return grid;
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

  const Eigen::MatrixXd& mo_coefficients =_orbitals.MOs().eigenvectors();
  const Eigen::MatrixXd& mo_energies = _orbitals.MOs().eigenvalues();
  const Eigen::MatrixXd overlap = OverlapMatrix();
  Eigen::MatrixXd H = overlap * mo_coefficients * mo_energies *
                      mo_coefficients.transpose() * overlap;

  return H;
}
Eigen::MatrixXcd SternheimerW::CoulombMatrix(Eigen::Vector3d gridpoint) const{

  AOBasis basis = _orbitals.SetupDftBasis();
  AOMultipole aoesp;
  aoesp.FillPotential(basis,gridpoint);
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

Eigen::MatrixXcd SternheimerW::DeltaNOneShot(
    std::complex<double> w, Eigen::Vector3d r) const{

  Eigen::MatrixXcd V = CoulombMatrix(r);
  
  double alpha = 2 * (_mo_energies(_num_occ_lvls) - _mo_energies(2));

  Eigen::MatrixXcd solution_p =
      Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);
  Eigen::MatrixXcd solution_m =
      Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);

  Eigen::MatrixXcd H_new;
  Eigen::MatrixXcd LHS_P;
  Eigen::MatrixXcd LHS_M;

  Eigen::VectorXcd RHS;// = Eigen::Xcd::Zero(_basis_size, _basis_size);

  //Multishift::MultiShiftResult result;
  
  for (int v = 0; v < _num_occ_lvls; v++) {

    RHS = SternheimerRHS(_inverse_overlap, _density_Matrix, V, _mo_coefficients.col(v));

      
        LHS_P = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, true);
        solution_p.col(v) = LHS_P.colPivHouseholderQr().solve(RHS);
        //solution_p[i].col(v) = _multishift.DoMultishift(LHS_P,RHS,w[i],result);
//        if(((LHS_P+w[i]*Eigen::MatrixXcd::Identity(_basis_size,_basis_size))*solution_p[i].col(v)-RHS).norm()>1e-13){
//            std::cout<<"res_p="<<(LHS_P+w[i]*Eigen::MatrixXcd::Identity(_basis_size,_basis_size))*solution_p[i].col(v)-RHS<<std::endl;
//        }

        LHS_M = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap, _mo_energies(v), w, false);
        
        solution_m.col(v) = LHS_M.colPivHouseholderQr().solve(RHS);
        //solution_m[i].col(v) = _multishift.DoMultishift(LHS_M,RHS,-w[i],result);
        //if(((LHS_M+w[i]*Eigen::MatrixXcd::Identity(_basis_size,_basis_size))*solution_m[i].col(v)-RHS).norm()>1e-13){
        //    std::cout<<"res_m="<<(LHS_M-w[i]*Eigen::MatrixXcd::Identity(_basis_size,_basis_size))*solution_m[i].col(v)-RHS<<std::endl;
        //}
      
    }
  

  Eigen::MatrixXcd delta_n=Eigen::MatrixXcd::Zero(solution_p.cols(),solution_p.rows());

    delta_n +=
    2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) * solution_p.transpose() +
    2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) * solution_m.transpose();

  return delta_n;
}

Eigen::MatrixXcd DielectricMatrix(Eigen::MatrixXcd deltaN){

  return Eigen::MatrixXcd::Identity(deltaN.rows(),deltaN.cols())-deltaN;

}

Eigen::MatrixXcd SternheimerW::GreensFunctionLHS(std::complex<double> w) const{

  return _Hamiltonian_Matrix-w*_overlap_Matrix;


}

Eigen::MatrixXcd SternheimerW::AnalyticGreensfunction(std::complex<double> w, Eigen::Vector3d r) const{

  return GreensFunctionLHS(w).colPivHouseholderQr().solve(-1*Eigen::MatrixXcd::Identity(_basis_size,_basis_size));

}

std::vector<Eigen::MatrixXcd> SternheimerW::Polarisability(
    double omega_start, double omega_end, int steps, double imaginary_shift,
    double lorentzian_broadening, int resolution_output, std::string gridtype) {
  
  std::vector<std::complex<double>> grid_w =
      BuildGrid(omega_start, omega_end, steps, imaginary_shift);
  std::vector<std::complex<double>> w =
      BuildGrid(omega_start, omega_end, steps, 0);
  
  //initializePade(3);
  //initializeMultishift(_orbitals.getBasisSetSize());

  AOBasis dftbasis = _orbitals.SetupDftBasis();
  NumericalIntegration numint;
  numint.GridSetup(gridtype, _orbitals.QMAtoms(),
                   dftbasis);  // For now use medium grid
  std::vector<std::pair<double, const Eigen::Vector3d*>> grid =
      numint.getGridpoints();

  std::cout << "gridsize= " << numint.getGridSize() << std::endl;

  AOBasis basis = _orbitals.SetupDftBasis();
  
  std::vector<Eigen::MatrixXcd> Polar;
  
  int index=0;
  for (const std::pair<double, const Eigen::Vector3d*> point: grid) {  
    //std::vector<Eigen::MatrixXcd> delta_n = DeltaNOneShot(grid_w, *point.second);
    Eigen::Vector3d gridcont = *point.second;
    for (int o = 0; o < grid_w.size(); o++) {
      Eigen::MatrixXcd dielectricMatrix = DielectricMatrix(DeltaNOneShot(grid_w[o], *point.second));
    }
    if(index%1000==0){
        std::cout<<"Iteration="<<index<<std::endl;
    }
    index++;
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

}  // namespace xtp
}  // namespace votca