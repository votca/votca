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
#include <votca/xtp/sternheimer.h>

namespace votca {
namespace xtp {

Eigen::MatrixXd Sternheimer::CalculateOverlapMatrix() {

  Orbitals orbitals = this->_orbitals;

  const int& num_occ_lvls = orbitals.getNumberOfAlphaElectrons();
  AOBasis basis = orbitals.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  return overlap.Matrix();
}

Eigen::MatrixXd Sternheimer::CalculateDensityMatrix() {
  Orbitals orbitals = this->_orbitals;
  return orbitals.DensityMatrixGroundState();
}

Eigen::MatrixXd Sternheimer::CalculateHamiltonian() {

  Orbitals orbitals = this->_orbitals;
  const Eigen::MatrixXd& mo_coefficients = orbitals.MOCoefficients();
  const Eigen::MatrixXd& mo_energies = orbitals.MOEnergies().asDiagonal();
  const Eigen::MatrixXd overlap = CalculateOverlapMatrix();
  Eigen::MatrixXd H = overlap * mo_coefficients * mo_energies *
                      mo_coefficients.transpose() * overlap;

  Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
  ges.compute(H, overlap);

  //        std::cout<<"Overlap"<<std::endl;
  //        std::cout<<overlap<<std::endl;
  //        std::cout<<"Test Hamiltonian"<<std::endl;
  //        std::cout<<"Eigenvalues Hamiltonian"<<std::endl;
  //        std::cout<<ges.eigenvalues()<<std::endl;
  //        std::cout<<"MO_Energies"<<std::endl;
  //        std::cout<<mo_energies<<std::endl;
  //        std::cout<<"1 EigV"<<std::endl;
  //        std::cout<<mo_coefficients.col(0)<<std::endl;
  //        std::cout<<"1"<<std::endl;
  //        std::cout<<ges.eigenvectors().col(0)<<std::endl;

  // Eigen::MatrixXd H =
  // mo_coefficients.block(0,0,mo_coefficients.rows(),mo_coefficients.rows()) *
  // mo_energies.block(0,0,mo_coefficients.rows(),mo_coefficients.rows()) *
  // mo_coefficients.block(0,0,mo_coefficients.rows(),mo_coefficients.rows()).transpose();
  return H;
}

Eigen::MatrixXd Sternheimer::CalculateCoulombMatrix() {

  Orbitals orbitals = this->_orbitals;
  AOBasis basis = orbitals.SetupDftBasis();
  AOCoulomb coulomb;
  coulomb.Fill(basis);
  return coulomb.Matrix();
}

Eigen::MatrixXcd Sternheimer::SternheimerLHS(Eigen::MatrixXcd hamiltonian,
                                             Eigen::MatrixXcd overlap,
                                             double eps,
                                             std::complex<double> omega,
                                             bool pm) {

  Eigen::MatrixXcd S;  //= Eigen::MatrixXd::Zero(overlap.cols(),
                       //overlap.rows());
  // distinguish between +w and -w
  if (pm == true) {
    S = (eps + omega) * overlap;
  } else {
    S = (eps - omega) * overlap;
  }
  return (hamiltonian - S);
}

Eigen::VectorXcd Sternheimer::SternheimerRHS(Eigen::MatrixXcd overlap,
                                             Eigen::MatrixXcd density,
                                             Eigen::MatrixXcd pertubation,
                                             Eigen::VectorXcd coeff) {
  // Setup Identity Matrix
  Eigen::MatrixXcd I =
      Eigen::MatrixXcd::Identity(overlap.rows(), overlap.cols());
  // Perform Matrix Operati

  Eigen::VectorXcd M = -1 * (I - overlap * density) * pertubation * coeff;

  return M;
}

Eigen::VectorXcd Sternheimer::SternheimerSolve(Eigen::MatrixXcd& LHS,
                                               Eigen::VectorXcd& RHS) {

  Eigen::VectorXcd x;
  x = LHS.colPivHouseholderQr().solve(RHS);

  //        std::cout<<"error= "<<LHS*x-RHS<<std::endl;

  return x;

  //         Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  //
  //         //Solving Sternheimer equation using Biconjugate Gradients Method
  //         solver.compute(LHS);
  //         return solver.solve(RHS);
}

Eigen::MatrixXcd Sternheimer::DeltaVMatrix(Eigen::MatrixXcd deltaV) {

  Orbitals orbitals = this->_orbitals;

  const Eigen::MatrixXd& s = CalculateOverlapMatrix();

  Eigen::MatrixXcd deltaVEx =
      Eigen::MatrixXcd::Zero(deltaV.cols(), deltaV.rows());

  for (int i = 0; i < deltaV.cols(); i++) {
    for (int j = 0; j < deltaV.rows(); j++) {
      for (int a = 0; a < deltaV.rows(); a++) {
        for (int b = 0; b < deltaV.rows(); b++) {
          deltaVEx(i, j) = deltaVEx(i, j) + s(i, a) * deltaV(a, b) * s(b, j);
        }
      }
    }
  }

  return deltaVEx;
}

Eigen::MatrixXcd Sternheimer::CalculateDeltaN(std::complex<double> w) {

  Orbitals orbitals = this->_orbitals;

  const int& num_occ_lvls = orbitals.getNumberOfAlphaElectrons();
  const int& basis_size = orbitals.getBasisSetSize();
  // Setting up Matrices needed for the Sternheimer equation

  Eigen::MatrixXd H = CalculateHamiltonian();

  Eigen::MatrixXd S = CalculateOverlapMatrix();
  Eigen::MatrixXd p = CalculateDensityMatrix();
  Eigen::MatrixXd V = CalculateCoulombMatrix();
  const Eigen::MatrixXd& mo_coefficients = orbitals.MOCoefficients();

  const Eigen::VectorXd& mo_energies = orbitals.MOEnergies();

  // Initialising Solution Vectors and density/dielectric  matrix

  Eigen::MatrixXcd delta_c_p = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::MatrixXcd delta_c_m = Eigen::MatrixXcd::Zero(basis_size, basis_size);

  Eigen::MatrixXcd delta_n = Eigen::MatrixXcd::Zero(basis_size, basis_size);

  // Setting up Sternheimer equation and solving it

  Eigen::MatrixXcd LHS_P = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::MatrixXcd LHS_M = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::VectorXcd RHS = Eigen::VectorXcd::Zero(basis_size);
  Eigen::MatrixXcd H_new = Eigen::MatrixXcd::Zero(basis_size, basis_size);

  int homolvl = orbitals.getHomo();
  double alpha = 2 * (mo_energies(num_occ_lvls) - mo_energies(2));

  for (int v = 0; v < num_occ_lvls; v++) {
    if (w == 0.0) {
      H_new = H + alpha * S * p.transpose();
      LHS_P = SternheimerLHS(H_new, S, mo_energies(v), w, true);
      LHS_M = SternheimerLHS(H_new, S, mo_energies(v), w, false);
    } else {
      LHS_P = SternheimerLHS(H, S, mo_energies(v), w, true);
      LHS_M = SternheimerLHS(H, S, mo_energies(v), w, false);
    }

    RHS = SternheimerRHS(S, p, V, mo_coefficients.col(v));

    delta_c_p.col(v) = SternheimerSolve(LHS_P, RHS);
    delta_c_m.col(v) = SternheimerSolve(LHS_M, RHS);

    // calculating perturbed Density matrix

    //            delta_n+=mo_coefficients.col(v).transpose()*delta_c_p.col(v);
    //            delta_n+=mo_coefficients.col(v).transpose()*delta_c_m.col(v);
    //
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

  std::complex<double> N = (delta_n * S).trace();
  std::complex<double> N2 = (p * S).trace();

  std::cout << "Total change: " << N << std::endl;
  std::cout << "Total charge: " << N2 << std::endl;

  return delta_n;
}

Eigen::MatrixXcd Sternheimer::CalculateDielectricMatrix(
    std::complex<double> w) {

  Eigen::MatrixXd S = CalculateOverlapMatrix();
  Eigen::MatrixXcd delta_n = CalculateDeltaN(w);

  // Calculating the dielectric Matrix
  return S - S * delta_n * S;
}

Eigen::MatrixXcd Sternheimer::CalculateScreenedCoulombOS(
    std::complex<double> w) {

  Orbitals orbitals = this->_orbitals;

  Eigen::MatrixXd S = CalculateOverlapMatrix();
  Eigen::MatrixXd v = CalculateCoulombMatrix();
  Eigen::MatrixXcd eps = CalculateDielectricMatrix(w);

  return eps.inverse() * S * v;
}

Eigen::MatrixXcd Sternheimer::CalculateScreenedCoulombSC(
    std::complex<double> w) {

  Orbitals orbitals = this->_orbitals;

  // Setting up constants
  const int& num_occ_lvls = orbitals.getNumberOfAlphaElectrons();
  const int& basis_size = orbitals.getBasisSetSize();

  // Setting up needed matrices
  Eigen::MatrixXd overlap = CalculateOverlapMatrix();
  Eigen::MatrixXd density = CalculateDensityMatrix();
  Eigen::MatrixXd coulomb = CalculateCoulombMatrix();
  Eigen::MatrixXd Hamiltonian = CalculateHamiltonian();
  const Eigen::MatrixXd& mo_coefficients = orbitals.MOCoefficients();
  const Eigen::VectorXd& mo_energies = orbitals.MOEnergies();

  // Performing first loop with initial guess
  Eigen::MatrixXcd delta_n = CalculateDeltaN(w);

  Eigen::MatrixXcd delta_V = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  for (int i = 0; i < basis_size; i++) {
    for (int j = 0; j < basis_size; j++) {
      for (int n = 0; n < basis_size; n++) {
        for (int m = 0; m < basis_size; m++) {
          delta_V(i, m) = delta_n(i, j) * overlap(j, n) * coulomb(n, m);
        }
      }
    }
  }

  Eigen::MatrixXcd W_new = coulomb + delta_V;

  // Initialising Matrices
  Eigen::MatrixXcd LHS_P = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::MatrixXcd LHS_M = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::VectorXcd RHS = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::MatrixXcd W_old = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::MatrixXcd H_new = Eigen::MatrixXcd::Zero(basis_size, basis_size);

  Eigen::MatrixXcd delta_c_p = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::MatrixXcd delta_c_m = Eigen::MatrixXcd::Zero(basis_size, basis_size);

  // Setting up max interations and tolerance
  double diff = 10000;
  int count = 0;
  double tol = 0.0001;
  int max_iter = 1000;

  double alpha = 2 * (mo_energies(num_occ_lvls) - mo_energies(0));

  while (diff > tol && count < max_iter) {

    // storing new W
    W_old = W_new;

    for (int v = 0; v < num_occ_lvls; v++) {

      if (w == 0.0) {

        H_new = Hamiltonian + alpha * overlap * density.transpose();
        LHS_P = SternheimerLHS(H_new, overlap, mo_energies(v), w, true);
        LHS_M = SternheimerLHS(H_new, overlap, mo_energies(v), w, false);
      } else {
        LHS_P = SternheimerLHS(Hamiltonian, overlap, mo_energies(v), w, true);
        LHS_M = SternheimerLHS(Hamiltonian, overlap, mo_energies(v), w, false);
      }

      RHS = SternheimerRHS(overlap, density, DeltaVMatrix(delta_V) + coulomb,
                           mo_coefficients.col(v));

      delta_c_p.col(v) = SternheimerSolve(LHS_P, RHS);
      delta_c_m.col(v) = SternheimerSolve(LHS_M, RHS);

      // calculating perturbed Density matrix
    }
    // calculating new W
    delta_n = Eigen::MatrixXd::Zero(basis_size, basis_size);
    delta_V = Eigen::MatrixXd::Zero(basis_size, basis_size);

    for (int v = 0; v < num_occ_lvls; v++) {
      for (int i = 0; i < basis_size; i++) {
        for (int j = 0; j < basis_size; j++) {
          delta_n(i, j) = delta_n(i, j) +
                          2 * mo_coefficients(i, v) * delta_c_p(j, v) +
                          2 * mo_coefficients(i, v) * delta_c_m(j, v);
        }
      }
    }
    for (int i = 0; i < basis_size; i++) {
      for (int j = 0; j < basis_size; j++) {
        for (int n = 0; n < basis_size; n++) {
          for (int m = 0; m < basis_size; m++) {
            delta_V(i, m) = delta_n(i, j) * overlap(j, n) * coulomb(n, m);
          }
        }
      }
    }

    //            std::cout<<"deltaV"<<std::endl;
    //            std::cout<<delta_V<<std::endl;

    W_new = coulomb + delta_V;

    count++;
    diff = (W_new - W_old).norm();

    std::cout << "iteration: " << count << std::endl;
    std::cout << "diff: " << diff << std::endl;

    if (diff < tol) {
      std::cout << "Self-consistency reached after " << count << "iterations."
                << std::endl;
      std::cout << "The last difference was " << diff << "." << std::endl;
    }
    if (count > max_iter - 1) {
      std::cout << "Max interations reached with " << count << "iterations."
                << std::endl;
      std::cout << "The last difference was " << diff << "." << std::endl;
    }

    delta_n = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  }

  return W_new;
}

Eigen::VectorXcd Sternheimer::PadeAppox(Eigen::VectorXcd grid) {

  Orbitals orbitals = this->_orbitals;
  const int& basis_size = orbitals.getBasisSetSize();

  Eigen::MatrixXcd W = Eigen::MatrixXcd::Zero(basis_size, basis_size);

  for (int i = 0; i < grid.size(); i++) {

    W = CalculateScreenedCoulombSC(grid(i));
  }
}

Eigen::MatrixXcd Sternheimer::NonanalyticGreensfunction(Orbitals& orbitals,
                                                        double w) {

  const int& basis_size = orbitals.getBasisSetSize();
  const int& num_occ_lvls = orbitals.getNumberOfAlphaElectrons();

  const Eigen::VectorXd& mo_energies = orbitals.MOEnergies();
  const Eigen::MatrixXd& mo_coefficients = orbitals.MOCoefficients();

  Eigen::MatrixXcd Gn = Eigen::MatrixXcd::Zero(basis_size, basis_size);
  Eigen::MatrixXcd add = Eigen::MatrixXcd::Zero(basis_size, basis_size);

  const double pi = 3.14159265358979323846;

  std::complex<double> c(0, 2 * pi);

  for (int v = 0; v < num_occ_lvls; v++) {

    if (mo_energies(v) == w) {

      add = c * mo_coefficients.col(v) * mo_coefficients.col(v).transpose();

      Gn = Gn + add;
    }
  }
  return Gn;
}

Eigen::MatrixXcd Sternheimer::AnalyticGreensLHS(Orbitals& orbitals, double w,
                                                double eta) {

  Eigen::MatrixXd H = CalculateHamiltonian();
  Eigen::MatrixXd S = CalculateOverlapMatrix();

  std::complex<double> c(w, eta);

  return H - c * S;
}

Eigen::MatrixXcd Sternheimer::AnalyticGreensRHS(Orbitals& orbitals, double w) {

  const int& num_occ_lvls = orbitals.getNumberOfAlphaElectrons();
  const int& basis_size = orbitals.getBasisSetSize();

  return Eigen::MatrixXd::Identity(basis_size, basis_size);
}

Eigen::MatrixXcd Sternheimer::AnalyticGreensolver(Eigen::MatrixXcd A,
                                                  Eigen::MatrixXcd b) {

  Eigen::MatrixXcd x;
  x = A.colPivHouseholderQr().solve(b);

  return x;
}

Eigen::MatrixXcd Sternheimer::AnalyticGreensfunction(Orbitals& orbitals,
                                                     double w, double eta) {

  Eigen::MatrixXcd LHS = AnalyticGreensLHS(orbitals, w, eta);
  Eigen::MatrixXcd RHS = AnalyticGreensRHS(orbitals, w);

  return AnalyticGreensolver(LHS, RHS);
}

Eigen::MatrixXcd Sternheimer::CalculateGreensfunction(Orbitals& orbitals,
                                                      double w, double eta) {

  Eigen::MatrixXcd G_A = AnalyticGreensfunction(orbitals, w, eta);
  Eigen::MatrixXcd G_N = NonanalyticGreensfunction(orbitals, w);

  return G_A + G_N;
}

Eigen::MatrixXcd Sternheimer::SelfEnergyEx(Orbitals orbitals) {

  const int& num_occ_lvls = orbitals.getNumberOfAlphaElectrons();

  Eigen::MatrixXd V = CalculateCoulombMatrix();
  const Eigen::MatrixXd& mo_coefficients = orbitals.MOCoefficients();

  Eigen::MatrixXd SE =
      Eigen::MatrixXd::Zero(mo_coefficients.rows(), mo_coefficients.cols());

  for (int v = 0; v < num_occ_lvls; v++) {

    SE += mo_coefficients.col(v).transpose() * mo_coefficients.col(v) * V;
  }

  return -SE;
}

Eigen::MatrixXcd Sternheimer::SelfEnergyC(Orbitals orbitals,
                                          Eigen::VectorXd grid, double w,
                                          double eta) {

  const double pi = 3.14159265358979323846;

  Eigen::MatrixXd V = CalculateCoulombMatrix();

  const Eigen::MatrixXd& mo_coefficients = orbitals.MOCoefficients();

  Eigen::MatrixXcd SE =
      Eigen::MatrixXd::Zero(mo_coefficients.rows(), mo_coefficients.cols());

  double step = 0.0;

  for (int i = 0; i < grid.size(); i++) {

    step = grid(i + 1) - grid(i);

    SE += AnalyticGreensfunction(orbitals, w + grid(i), eta) *
          (CalculateScreenedCoulombOS(grid(i)) - V) * step;
  }
  std::complex<double> c(0, 1 / (2 * pi));

  return c * SE;
}

Eigen::MatrixXcd Sternheimer::CalculateSelfEnergy(Orbitals& orbitals, double w,
                                                  double eta,
                                                  Eigen::VectorXd grid) {

  Eigen::MatrixXcd SEC = SelfEnergyC(orbitals, grid, w, eta);
  Eigen::MatrixXcd SEEx = SelfEnergyEx(orbitals);

  //        std::cout<<"SEC"<<std::endl;
  //        std::cout<<SEC<<std::endl;
  //        std::cout<<"SEEx"<<std::endl;
  //        std::cout<<SEEx<<std::endl;

  return SEC + SEEx;
}

//    double Sternheimer::CalculateSpectralfunction(Eigen::MatrixXcd Selfenergy,
//    Eigen::MatrixXcd deltaSE, double w){
//
//        const double pi = 3.14159265358979323846;
//
//        const Eigen::VectorXd& mo_energies = orbitals.MOEnergies();
//
//        double A = 0;
//
//        for(int n=0;n<Selfenergy.rows();n++){
//
//            A+=(std::abs(std::imag(Selfenergy(n,0))))/(pow(w-mo_energies(n)-std::real(deltaSE(n,0)),2)+pow(std::imag(Selfenergy(n,0)),2));
//
//        }
//
//        return (1/pi)*A;
//
//    }
//
//    Eigen::VectorXcd Sternheimer::CalculateqpCorrection(Eigen::MatrixXcd
//    Pertubation, Orbitals& orb){
//
//        const Eigen::MatrixXd& mo_coefficients = orb.MOCoefficients();
//
//        Eigen::VectorXcd qpc = Eigen::VectorXcd::Zero(mo_coefficients.cols());
//
//        for(int n=0;n<mo_coefficients.cols();n++){
//
//            qpc(n)=mo_coefficients.col(n).transpose()*Pertubation*mo_coefficients.col(n);
//
//        }
//        return qpc;
//    }
//
//
//    Eigen::VectorXcd Sternheimer::CalculateqpCorrection(Eigen::MatrixXcd
//    SelfEnergy, Eigen::MatrixXcd VXC, Orbitals& orb){
//
//        const Eigen::MatrixXd& mo_coefficients = orb.MOCoefficients();
//
//        Eigen::VectorXcd qpc = Eigen::VectorXcd::Zero(mo_coefficients.cols());
//
//        for(int n=0;n<mo_coefficients.cols();n++){
//
//            qpc(n)=mo_coefficients.col(n).transpose()*(SelfEnergy-VXC)*mo_coefficients.col(n);
//
//
//        }
//        return qpc;
//
//    }

Eigen::MatrixXcd Sternheimer::CalculateKSselfenergie(
    Eigen::MatrixXcd Selfenergie, Eigen::MatrixXd XC) {}

Eigen::MatrixXd DipoleMatrix() {}

bool Sternheimer::evaluate() { return true; }

//    Eigen::MatrixXd Sternheimer::Polarisability(Orbitals orb,Eigen::MatrixXd
//    delta_n, int x, int y, double w){
//
//        Eigen::MatrixXd chi = CalculateDeltaN(orb, w);
//
//
//    }

}  // namespace xtp
}  // namespace votca
