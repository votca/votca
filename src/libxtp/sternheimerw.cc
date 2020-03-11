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
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aomatrix3d.h>
#include <votca/xtp/aopotential.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/multishift.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/padeapprox.h>
#include <votca/xtp/sternheimerw.h>

namespace votca {
namespace xtp {

Eigen::VectorXd SternheimerW::EvaluateBasisAtPosition(
    const AOBasis& dftbasis, const Eigen::Vector3d& pos) const {

  // get value of orbitals at each gridpoint
  Eigen::VectorXd tmat = Eigen::VectorXd::Zero(dftbasis.AOBasisSize());
  for (const AOShell& shell : dftbasis) {
    const double decay = shell.getMinDecay();
    const Eigen::Vector3d& shellpos = shell.getPos();
    Eigen::Vector3d dist = shellpos - pos;
    double distsq = dist.squaredNorm();
    // if contribution is smaller than -ln(1e-10), calc density
    if ((decay * distsq) < 20.7) {
      Eigen::VectorBlock<Eigen::VectorXd> tmat_block =
          tmat.segment(shell.getStartIndex(), shell.getNumFunc());
      shell.EvalAOspace(tmat_block, pos);
    }
  }
  return tmat;
}

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


void SternheimerW::configurate(const options_sternheimer& opt) { _opt = opt; }

std::vector<std::complex<double>> SternheimerW::BuildGrid() const {

  std::vector<std::complex<double>> grid;

  double stepsize = (_opt.start_frequency_grid - _opt.end_frequency_grid) / (_opt.number_of_frequency_grid_points);
  const double ev2hrt = votca::tools::conv::ev2hrt;
  std::complex<double> d(1, 0);
  std::complex<double> i(0, 1);

  for (Index n = 0; n < _opt.number_of_frequency_grid_points; n++) {
    // Converts from input eV to Hartree
    grid.push_back(_opt.start_frequency_grid * ev2hrt + n * stepsize * ev2hrt +
                  _opt.imaginary_shift * i * ev2hrt);
  }

  return grid;
}

void SternheimerW::initializeMultishift(Index size) {
  _multishift.setMatrixSize(size);
}
void SternheimerW::initializePade(Index size) { _pade.initialize(size); }

Eigen::MatrixXd SternheimerW::OverlapMatrix() const {

  AOBasis basis = _orbitals.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  return overlap.Matrix();
}

Eigen::MatrixXd SternheimerW::DensityMatrix() const {

  return _orbitals.DensityMatrixGroundState();
}

Eigen::MatrixXd SternheimerW::Hamiltonian() const {

  const Eigen::MatrixXd& mo_coefficients = _orbitals.MOs().eigenvectors();
  const Eigen::MatrixXd& mo_energies = _orbitals.MOs().eigenvalues();
  const Eigen::MatrixXd overlap = OverlapMatrix();
  Eigen::MatrixXd H = overlap * mo_coefficients * mo_energies *
                      mo_coefficients.transpose() * overlap;

  return H;
}
Eigen::MatrixXcd SternheimerW::CoulombMatrix(Eigen::Vector3d gridpoint) const {

  AOBasis basis = _orbitals.SetupDftBasis();
  AOMultipole aoesp;
  aoesp.FillPotential(basis, gridpoint);
  return aoesp.Matrix();
}

Eigen::MatrixXcd SternheimerW::SternheimerLHS(
    const Eigen::MatrixXcd& hamiltonian,
    const Eigen::MatrixXcd& inverse_overlap, double eps,
    std::complex<double> omega, bool pm) const {

  // distinguish between +w and -w
  std::complex<double> temp = eps + omega;
  if (pm != true) {
    temp = (eps - omega);
  }
  return (inverse_overlap * hamiltonian -
          temp * Eigen::MatrixXcd::Identity(_basis_size, _basis_size));
}

Eigen::VectorXcd SternheimerW::SternheimerRHS(
    const Eigen::MatrixXcd& inverse_overlap, const Eigen::MatrixXcd& density,
    const Eigen::MatrixXcd& pertubation, const Eigen::VectorXcd& coeff) const {

  return -1 * (inverse_overlap - density) * pertubation * coeff;
}

Eigen::VectorXcd SternheimerW::SternheimerSolve(const Eigen::MatrixXcd& LHS,
                                                const Eigen::VectorXcd& RHS) {
  return LHS.colPivHouseholderQr().solve(RHS);
}

Eigen::MatrixXcd SternheimerW::DeltaNOneShot(std::complex<double> w,
                                             Eigen::Vector3d r) const {

  Eigen::MatrixXcd V = CoulombMatrix(r);

  Eigen::MatrixXcd solution_p =
      Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);
  Eigen::MatrixXcd solution_m =
      Eigen::MatrixXcd::Zero(_basis_size, _num_occ_lvls);

  Eigen::MatrixXcd H_new;
  Eigen::MatrixXcd LHS_P;
  Eigen::MatrixXcd LHS_M;

  Eigen::VectorXcd RHS;  // = Eigen::Xcd::Zero(_basis_size, _basis_size);

  for (Index v = 0; v < _num_occ_lvls; v++) {

    RHS = SternheimerRHS(_inverse_overlap, _density_Matrix, V,
                         _mo_coefficients.col(v));
    LHS_P = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap,
                           _mo_energies(v), w, true);
    solution_p.col(v) = LHS_P.colPivHouseholderQr().solve(RHS);
    LHS_M = SternheimerLHS(_Hamiltonian_Matrix, _inverse_overlap,
                           _mo_energies(v), w, false);
    solution_m.col(v) = LHS_M.colPivHouseholderQr().solve(RHS);
  }

  Eigen::MatrixXcd delta_n =
      2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
          solution_p.transpose() +
      2 * _mo_coefficients.block(0, 0, _basis_size, _num_occ_lvls) *
          solution_m.transpose();

  return delta_n;
}

Eigen::MatrixXcd SternheimerW::DielectricMatrix(Eigen::MatrixXcd deltaN) const {
  deltaN.diagonal().array() -= 1.0;
  return -deltaN;
}
void SternheimerW::printDielectricFunction(double omega_start, double omega_end,
                                           Index steps, double imaginary_shift,
                                           double lorentzian_broadening,
                                           Index resolution_output,
                                           std::string gridtype) {

  std::vector<std::complex<double>> grid_w =
      BuildGrid();

  AOBasis dftbasis = _orbitals.SetupDftBasis();

  // NumericalIntegration numint;
  // numint.GridSetup(gridtype, _orbitals.QMAtoms(), dftbasis);
  std::vector<std::pair<double, const Eigen::Vector3d*>> grid;

  for (Index o = 0; o < grid_w.size(); o++) {
    for (const std::pair<double, const Eigen::Vector3d*> point : grid) {
      Eigen::Vector3d gridcont = *point.second;
      Eigen::MatrixXcd dielectricMatrix =
          DielectricMatrix(DeltaNOneShot(grid_w[o], *point.second));
      Eigen::MatrixXcd ScreenedCoulomb =
          dielectricMatrix.inverse() * CoulombMatrix(*point.second);
    }
  }
}

std::complex<double> SternheimerW::ScreenedCoulomb(
    Eigen::Vector3d gridpoint1, Eigen::Vector3d gridpoint2,
    std::complex<double> frequency) const {

  AOBasis basis = _orbitals.SetupDftBasis();
  Eigen::MatrixXcd dielectricMatrix =
      DielectricMatrix(DeltaNOneShot(frequency, gridpoint1));
  Eigen::MatrixXcd ScreenedCoulombMatrix =
      dielectricMatrix.inverse() * CoulombMatrix(gridpoint1);
  std::complex<double> W = 0.0;

  for (Index n = 0; n < _basis_size; n++) {
    for (Index m = 0; m < _basis_size; m++) {

      W += EvaluateBasisAtPosition(basis, gridpoint2)(n) *
           ScreenedCoulombMatrix(n, m) *
           EvaluateBasisAtPosition(basis, gridpoint2)(m);
    }
  }
  return W;
}

Eigen::MatrixXcd SternheimerW::ScreenedCoulomb(
    Eigen::Vector3d gridpoint1, std::complex<double> frequency) const {
  AOBasis basis = _orbitals.SetupDftBasis();
  Eigen::MatrixXcd dielectricMatrix =
      DielectricMatrix(DeltaNOneShot(frequency, gridpoint1)).inverse();
  dielectricMatrix.diagonal().array() -= 1.0; // We only want the correlation part
     
   Eigen::MatrixXcd screened =  dielectricMatrix * CoulombMatrix(gridpoint1);
   
  return screened;
  
}

Eigen::MatrixXcd SternheimerW::GreensFunctionLHS(std::complex<double> w) const {

  return _Hamiltonian_Matrix - w * _overlap_Matrix;
}

Eigen::MatrixXcd SternheimerW::AnalyticGreensfunction(
    std::complex<double> w) const {
  return GreensFunctionLHS(w).colPivHouseholderQr().solve(
      -1 * Eigen::MatrixXcd::Identity(_basis_size, _basis_size));
}

double SternheimerW::Lorentzian(double center,
                                std::complex<double> freq) const {
  double gamma = 0.001;

  // double prefactor = 1/(votca::tools::conv::Pi*gamma);
  double enumerator = gamma * gamma;
  double denominator = (std::pow(abs(freq - center), 2) + gamma * gamma);
  return (enumerator / denominator);
}

Eigen::MatrixXcd SternheimerW::NonAnalyticGreensfunction(
    std::complex<double> freq) const {
  Eigen::MatrixXcd NonAnalyticGreens =
      Eigen::MatrixXcd::Zero(_basis_size, _basis_size);

  std::complex<double> factor(0, 2 * votca::tools::conv::Pi);

  for (Index v = 0; v < _num_occ_lvls; v++) {
    NonAnalyticGreens += factor * Lorentzian(_mo_energies[v], freq) *
                         _mo_coefficients.col(v) *
                         _mo_coefficients.col(v).transpose();
  }

  return NonAnalyticGreens;
}

std::complex<double> SternheimerW::GreensFunction(
    Eigen::Vector3d gridpoint1, Eigen::Vector3d gridpoint2,
    std::complex<double> frequency) const {

  AOBasis basis = _orbitals.SetupDftBasis();

  Eigen::MatrixXcd GreensFunctionMatrix =
      NonAnalyticGreensfunction(frequency) + AnalyticGreensfunction(frequency);

  std::complex<double> G = 0.0;

  for (Index n = 0; n < _basis_size; n++) {
    for (Index m = 0; m < _basis_size; m++) {

      G += EvaluateBasisAtPosition(basis, gridpoint1)(n) *
           GreensFunctionMatrix(n, m) *
           EvaluateBasisAtPosition(basis, gridpoint2)(m);
    }
  }
  return G;
}

Eigen::MatrixXcd SternheimerW::GreensFunction(
    std::complex<double> frequency) const {

  Eigen::MatrixXcd GreensFunctionMatrix =
      NonAnalyticGreensfunction(frequency) + AnalyticGreensfunction(frequency);
  return GreensFunctionMatrix;
}

void SternheimerW::printGreensfunction(double omega_start, double omega_end,
                                       Index steps, double imaginary_shift,
                                       double lorentzian_broadening,
                                       Index resolution_output) {

  std::vector<std::complex<double>> grid_w =
      BuildGrid();

  for (std::complex<double> w : grid_w) {

    double Greens = imag(
        (AnalyticGreensfunction(w) + NonAnalyticGreensfunction(w)).trace());

    std::cout << real(w) << " " << Greens << " "
              << real(AnalyticGreensfunction(w).trace()) << " "
              << imag(NonAnalyticGreensfunction(w).trace()) << std::endl;
  }
}

std::complex<double> SternheimerW::SelfEnergy_at_r(
    double omega, Eigen::Vector3d gridpoint1, Index n, Index np) const {

  std::vector<std::complex<double>> w =
      BuildGrid();

  AOBasis basis = _orbitals.SetupDftBasis();
  Eigen::VectorXd chi = EvaluateBasisAtPosition(basis, gridpoint1);
  Eigen::MatrixXd left = chi * chi.transpose();
  Eigen::MatrixXcd right = Eigen::MatrixXcd::Zero(_basis_size, _basis_size);
  // Integral over frequencies omega: trapezoidal rule
  double delta = (_opt.end_frequency_grid - _opt.start_frequency_grid) / _opt.number_of_frequency_grid_points;
  Index i = 0;
  for (std::complex<double> omega_p : w) {
    std::complex<double> exponent (-(1e-3),0.);
    exponent *= omega_p;
    
    if ( i == 0 ){
      right += 0.5* GreensFunction(omega + omega_p) * ScreenedCoulomb(gridpoint1, omega_p) * std::exp(exponent);
    }
    else if ( i == _opt.number_of_frequency_grid_points -1 ){
      right += 0.5* GreensFunction(omega + omega_p) * ScreenedCoulomb(gridpoint1, omega_p) * std::exp(exponent);
    }
    else {
    right +=
        GreensFunction(omega + omega_p) * ScreenedCoulomb(gridpoint1, omega_p) * std::exp(exponent);
  }
  i++;
  }
  
  return _mo_coefficients.col(n).transpose() * left * right *
         _mo_coefficients.col(np);
}


}  // namespace xtp
}  // namespace votca