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

Eigen::VectorXd DIIS::CalcCoeff() {
  success = true;
  const Index size = Index(errormatrixhist_.size());

  // C2-DIIS
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(size, size);

  for (Index i = 0; i < B.rows(); i++) {
    for (Index j = 0; j <= i; j++) {
      B(i, j) = Diis_Bs_[i][j];
      if (i != j) {
        B(j, i) = B(i, j);
      }
    }
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
  Eigen::MatrixXd eigenvectors = Eigen::MatrixXd::Zero(size, size);

  for (Index i = 0; i < es.eigenvectors().cols(); i++) {
    double norm = es.eigenvectors().col(i).sum();
    eigenvectors.col(i) = es.eigenvectors().col(i) / norm;
  }

  // Choose solution by picking out solution with smallest absolute error
  Eigen::VectorXd errors =
      (eigenvectors.transpose() * B * eigenvectors).diagonal().cwiseAbs();

  double MaxWeight = 10.0;
  Index mincoeff = 0;
  success = false;
  for (Index i = 0; i < errors.size(); i++) {
    errors.minCoeff(&mincoeff);
    if (std::abs(eigenvectors.col(mincoeff).maxCoeff()) > MaxWeight) {
      errors[mincoeff] = std::numeric_limits<double>::max();
    } else {
      success = true;
      break;
    }
  }

  Eigen::VectorXd coeffs = eigenvectors.col(mincoeff);

  if (std::abs(coeffs[coeffs.size() - 1]) < 0.001) {
    success = false;
  }

  return coeffs;
}

}  // namespace xtp
}  // namespace votca
