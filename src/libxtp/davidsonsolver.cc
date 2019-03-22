/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <iostream>
#include <stdexcept>

#include <votca/xtp/davidsonsolver.h>
#include <votca/xtp/eigen.h>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

using namespace std;

DavidsonSolver::DavidsonSolver(ctp::Logger &log) : _log(log) {}

void DavidsonSolver::set_correction(std::string method) {
  if (method == "DPR")
    this->davidson_correction = CORR::DPR;
  else if (method == "OLSEN")
    this->davidson_correction = CORR::OLSEN;
  else
    throw std::runtime_error(method +
                             " is not a valid Davidson correction method");
}

Eigen::ArrayXi DavidsonSolver::_sort_index(Eigen::VectorXd &V) const {
  /* return the index of the sorted vector */
  Eigen::ArrayXi idx = Eigen::ArrayXi::LinSpaced(V.rows(), 0, V.rows() - 1);
  std::sort(idx.data(), idx.data() + idx.size(),
            [&](int i1, int i2) { return V[i1] < V[i2]; });
  return idx;
}

Eigen::MatrixXd DavidsonSolver::_get_initial_eigenvectors(
    Eigen::VectorXd &d, int size_initial_guess) const {

  /* Initialize the guess eigenvector so that they 'target' the lowest diagonal
   * elements */

  Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(d.size(), size_initial_guess);
  Eigen::ArrayXi idx = DavidsonSolver::_sort_index(d);

  for (int j = 0; j < size_initial_guess; j++) {
    guess(idx(j), j) = 1.0;
  }

  return guess;
}

Eigen::VectorXd DavidsonSolver::_dpr_correction(Eigen::VectorXd &r,
                                                Eigen::VectorXd &D,
                                                double lambda) const {
  /* Compute the diagonal preconditoned residue : delta = - (D - lambda)^{-1} r
   */
  
  Eigen::VectorXd delta = r.array() / (lambda - D.array());
  return delta;
}

Eigen::VectorXd DavidsonSolver::_olsen_correction(Eigen::VectorXd &r,
                                                  Eigen::VectorXd &x,
                                                  Eigen::VectorXd &D,
                                                  double lambda) const {
  /* Compute the olsen correction :

  \delta = (D-\lambda)^{-1} (-r + \epsilon x)

  */

  int size = r.rows();
  Eigen::VectorXd delta = Eigen::VectorXd::Zero(size);

  delta = DavidsonSolver::_dpr_correction(r, D, lambda);

  double _num = -x.transpose() * delta;
  double _denom =
      -x.transpose() * DavidsonSolver::_dpr_correction(x, D, lambda);
  double eps = _num / _denom;

  delta += eps * x;

  return delta;
}

Eigen::MatrixXd DavidsonSolver::_QR(Eigen::MatrixXd &A) const {

  int nrows = A.rows();
  int ncols = A.cols();
  ncols = std::min(nrows, ncols);

  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  Eigen::MatrixXd result =
      qr.householderQ() * Eigen::MatrixXd::Identity(nrows, ncols);
  return result;
}

Eigen::MatrixXd DavidsonSolver::_gramschmidt( Eigen::MatrixXd &A, int nstart ) const
{
    Eigen::MatrixXd Q = A;

    for(unsigned int j = nstart; j < A.cols(); ++j) {
        // Replace inner loop over each previous vector in Q with fast matrix-vector multiplication
        Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));
        // Normalize vector if possible (othw. means colums of A almsost lin. dep.
        if( Q.col(j).norm() <= 10e-14 * A.col(j).norm() ) {
            std::cerr << "Gram-Schmidt failed because A has lin. dep columns. Bye." << std::endl;
            break;
        } else {
            Q.col(j).normalize();
        }
    }
    return Q;
}

}  // namespace xtp
}  // namespace votca