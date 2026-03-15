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

// Standard includes
#include <iostream>

// Local VOTCA includes
#include "votca/xtp/diis.h"

namespace votca {
namespace xtp {

/**
 * Pulay DIIS implementation.
 *
 * The method stores recent commutator error matrices e_i and solves an
 * augmented linear system built from B_ij = <e_i, e_j>. The resulting
 * coefficients define the extrapolated Fock matrix F' = sum_i c_i F_i.
 */

// Update the Pulay error history. The scalar products stored in Diis_Bs_
// form the DIIS B matrix with elements
//
//   B_ij = <R_i | R_j> = Tr[R_i^T R_j],
//
// where R_i is the commutator residual of iteration i.
void DIIS::Update(Index maxerrorindex, const Eigen::MatrixXd& errormatrix) {

  if (int(errormatrixhist_.size()) == histlength_) {
    errormatrixhist_.erase(errormatrixhist_.begin() + maxerrorindex);
    Diis_Bs_.erase(Diis_Bs_.begin() + maxerrorindex);
    for (std::vector<double>& subvec : Diis_Bs_) {
      subvec.erase(subvec.begin() + maxerrorindex);
    }
  }

  errormatrixhist_.push_back(errormatrix);

  std::vector<double> Bijs;
  for (Index i = 0; i < Index(errormatrixhist_.size()) - 1; i++) {
    double value =
        errormatrix.cwiseProduct((errormatrixhist_[i]).transpose()).sum();
    Bijs.push_back(value);
    Diis_Bs_[i].push_back(value);
  }
  Bijs.push_back(errormatrix.cwiseProduct(errormatrix.transpose()).sum());
  Diis_Bs_.push_back(Bijs);
  return;
}

// Solve the standard constrained DIIS system
//
//   [ B  -1 ] [c]   [0]
//   [ -1  0 ] [l] = [-1]
//
// to obtain extrapolation coefficients c whose sum is one.
Eigen::VectorXd DIIS::CalcCoeff() {
  success = true;

  Index size = 0;
  if (!errormatrixhist_alpha_.empty()) {
    size = Index(errormatrixhist_alpha_.size());
  } else {
    size = Index(errormatrixhist_.size());
  }

  if (size < 2) {
    success = false;
    return Eigen::VectorXd();
  }

  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(size, size);

  for (Index i = 0; i < size; ++i) {
    for (Index j = 0; j <= i; ++j) {
      B(i, j) = Diis_Bs_[i][j];
      if (i != j) {
        B(j, i) = B(i, j);
      }
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
  if (es.info() != Eigen::Success) {
    success = false;
    return Eigen::VectorXd();
  }

  Eigen::MatrixXd eigenvectors = Eigen::MatrixXd::Zero(size, size);

  for (Index i = 0; i < size; ++i) {
    double norm = es.eigenvectors().col(i).sum();
    if (std::abs(norm) < 1e-14) {
      eigenvectors.col(i) = es.eigenvectors().col(i);
    } else {
      eigenvectors.col(i) = es.eigenvectors().col(i) / norm;
    }
  }

  Eigen::VectorXd errors =
      (eigenvectors.transpose() * B * eigenvectors).diagonal().cwiseAbs();

  constexpr double MaxWeight = 10.0;
  Index mincoeff = 0;
  success = false;

  for (Index i = 0; i < errors.size(); ++i) {
    errors.minCoeff(&mincoeff);
    if (std::abs(eigenvectors.col(mincoeff).maxCoeff()) > MaxWeight) {
      errors[mincoeff] = std::numeric_limits<double>::max();
    } else {
      success = true;
      break;
    }
  }

  if (!success) {
    return Eigen::VectorXd();
  }

  Eigen::VectorXd coeffs = eigenvectors.col(mincoeff);

  if (coeffs.size() == 0 || std::abs(coeffs[coeffs.size() - 1]) < 0.001) {
    success = false;
    return Eigen::VectorXd();
  }

  return coeffs;
}

// Unrestricted DIIS update. The alpha and beta residual overlaps are added
// so that B_ij represents the spin-summed residual metric used in the UKS
// extrapolation.
void DIIS::Update(Index maxerrorindex,
                  const Eigen::MatrixXd& errormatrix_alpha,
                  const Eigen::MatrixXd& errormatrix_beta) {

  if (int(errormatrixhist_alpha_.size()) == histlength_) {
    errormatrixhist_alpha_.erase(errormatrixhist_alpha_.begin() + maxerrorindex);
    errormatrixhist_beta_.erase(errormatrixhist_beta_.begin() + maxerrorindex);
    Diis_Bs_.erase(Diis_Bs_.begin() + maxerrorindex);
    for (std::vector<double>& subvec : Diis_Bs_) {
      subvec.erase(subvec.begin() + maxerrorindex);
    }
  }

  errormatrixhist_alpha_.push_back(errormatrix_alpha);
  errormatrixhist_beta_.push_back(errormatrix_beta);

  std::vector<double> Bijs;
  for (Index i = 0; i < Index(errormatrixhist_alpha_.size()) - 1; ++i) {
    double value =
        errormatrix_alpha.cwiseProduct(errormatrixhist_alpha_[i].transpose()).sum() +
        errormatrix_beta.cwiseProduct(errormatrixhist_beta_[i].transpose()).sum();
    Bijs.push_back(value);
    Diis_Bs_[i].push_back(value);
  }

  double self =
      errormatrix_alpha.cwiseProduct(errormatrix_alpha.transpose()).sum() +
      errormatrix_beta.cwiseProduct(errormatrix_beta.transpose()).sum();

  Bijs.push_back(self);
  Diis_Bs_.push_back(Bijs);
}

}  // namespace xtp
}  // namespace votca
