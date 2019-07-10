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
#include "votca/xtp/diis.h"
#include <iostream>

namespace votca {
namespace xtp {

void DIIS::Update(int maxerrorindex, const Eigen::MatrixXd& errormatrix) {

  if (int(_errormatrixhist.size()) == _histlength) {
    _errormatrixhist.erase(_errormatrixhist.begin() + maxerrorindex);
    _Diis_Bs.erase(_Diis_Bs.begin() + maxerrorindex);
    for (std::vector<double>& subvec : _Diis_Bs) {
      subvec.erase(subvec.begin() + maxerrorindex);
    }
  }

  _errormatrixhist.push_back(errormatrix);

  std::vector<double> Bijs;
  for (int i = 0; i < int(_errormatrixhist.size()) - 1; i++) {
    double value =
        errormatrix.cwiseProduct((_errormatrixhist[i]).transpose()).sum();
    Bijs.push_back(value);
    _Diis_Bs[i].push_back(value);
  }
  Bijs.push_back(errormatrix.cwiseProduct(errormatrix.transpose()).sum());
  _Diis_Bs.push_back(Bijs);
  return;
}

Eigen::VectorXd DIIS::CalcCoeff() {
  success = true;
  const bool _useold = false;
  const int size = _errormatrixhist.size();
  // old Pulat DIIs

  Eigen::VectorXd coeffs;
  if (_useold) {

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(size + 1, size + 1);
    Eigen::VectorXd a = Eigen::VectorXd::Zero(size + 1);
    a(size) = -1;
    for (int i = 0; i < B.rows() - 1; i++) {
      B(i, size) = -1;
      B(size, i) = -1;
    }

    for (int i = 0; i < B.rows() - 1; i++) {
      for (int j = 0; j <= i; j++) {
        B(i, j) = _Diis_Bs[i][j];
        if (i != j) {
          B(j, i) = B(i, j);
        }
      }
    }

    Eigen::VectorXd result = B.colPivHouseholderQr().solve(a);
    coeffs = result.head(size);
  } else {

    // C2-DIIS

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(size, size);

    for (int i = 0; i < B.rows(); i++) {
      for (int j = 0; j <= i; j++) {
        B(i, j) = _Diis_Bs[i][j];
        if (i != j) {
          B(j, i) = B(i, j);
        }
      }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(B);
    Eigen::MatrixXd eigenvectors = Eigen::MatrixXd::Zero(size, size);

    for (int i = 0; i < es.eigenvectors().cols(); i++) {
      double norm = es.eigenvectors().col(i).sum();
      eigenvectors.col(i) = es.eigenvectors().col(i) / norm;
    }

    // Choose solution by picking out solution with smallest absolute error
    Eigen::VectorXd errors =
        (eigenvectors.transpose() * B * eigenvectors).diagonal().cwiseAbs();

    double MaxWeight = 10.0;
    int mincoeff = 0;
    success = false;
    for (int i = 0; i < errors.size(); i++) {
      errors.minCoeff(&mincoeff);
      if (std::abs(eigenvectors.col(mincoeff).maxCoeff()) > MaxWeight) {
        errors[mincoeff] = std::numeric_limits<double>::max();
      } else {
        success = true;
        break;
      }
    }

    coeffs = eigenvectors.col(mincoeff);
  }

  if (std::abs(coeffs[coeffs.size() - 1]) < 0.001) {
    success = false;
  }

  return coeffs;
}

}  // namespace xtp
}  // namespace votca
