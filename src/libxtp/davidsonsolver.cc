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

void DavidsonSolver::set_ortho(std::string method) {
  if (method == "GS")
    this->_davidson_ortho = ORTHO::GS;
  else if (method == "QR")
    this->_davidson_ortho = ORTHO::QR;
  else
    throw std::runtime_error(
        method + " is not a valid Davidson orthogonalization method");
}

void DavidsonSolver::set_correction(std::string method) {
  if (method == "DPR")
    this->_davidson_correction = CORR::DPR;
  else if (method == "OLSEN")
    this->_davidson_correction = CORR::OLSEN;
  else
    throw std::runtime_error(method +
                             " is not a valid Davidson correction method");
}

void DavidsonSolver::set_tolerance(std::string tol) {
  if (tol == "loose")
    this->_tol = 1E-3;
  else if (tol == "normal")
    this->_tol = 1E-4;
  else if (tol == "strict")
    this->_tol = 1E-5;
  else
    throw std::runtime_error(tol + " is not a valid Davidson tolerance");
}

void DavidsonSolver::set_size_update(std::string update_size) {

  if (update_size == "min")
    this->_davidson_update = UPDATE::MIN;

  else if (update_size == "safe")
    this->_davidson_update = UPDATE::SAFE;

  else if (update_size == "max")
    this->_davidson_update = UPDATE::MAX;

  else
    throw std::runtime_error(update_size + " is not a valid Davidson update");
}

int DavidsonSolver::get_size_update(int neigen) {
  int size_update;
  switch (this->_davidson_update) {
    case UPDATE::MIN:
      size_update = neigen;
      break;
    case UPDATE::SAFE:
      if (neigen < 20)
        size_update = static_cast<int>(1.5 * neigen);
      else
        size_update = neigen + 10;
      break;
    case UPDATE::MAX:
      size_update = 2 * neigen;
      break;
  }
  return size_update;
}

Eigen::ArrayXi DavidsonSolver::argsort(Eigen::VectorXd &V) const {
  /* \brief return the index of the sorted vector */
  Eigen::ArrayXi idx = Eigen::ArrayXi::LinSpaced(V.rows(), 0, V.rows() - 1);
  std::sort(idx.data(), idx.data() + idx.size(),
            [&](int i1, int i2) { return V[i1] < V[i2]; });
  return idx;
}

Eigen::MatrixXd DavidsonSolver::SetupInitialEigenvectors(
    Eigen::VectorXd &d, int size_initial_guess) const {

  /* \brief Initialize the guess eigenvector so that they 'target' the lowest
   * diagonal elements */
  Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(d.size(), size_initial_guess);
  Eigen::ArrayXi idx = DavidsonSolver::argsort(d);

  for (int j = 0; j < size_initial_guess; j++) {
    guess(idx(j), j) = 1.0;
  }

  return guess;
}

Eigen::VectorXd DavidsonSolver::dpr_correction(Eigen::VectorXd &r,
                                               Eigen::VectorXd &D,
                                               double lambda) const {
  /* \brief Compute the diagonal preconditoned residue : delta = - (D -
   * lambda)^{-1} r
   */
  Eigen::VectorXd delta = r.array() / (lambda - D.array());
  return delta;
}

Eigen::VectorXd DavidsonSolver::olsen_correction(Eigen::VectorXd &r,
                                                 Eigen::VectorXd &x,
                                                 Eigen::VectorXd &D,
                                                 double lambda) const {
  /* \brief Compute the olsen correction :

  \delta = (D-\lambda)^{-1} (-r + \epsilon x)

  */
  int size = r.rows();
  Eigen::VectorXd delta = Eigen::VectorXd::Zero(size);
  delta = DavidsonSolver::dpr_correction(r, D, lambda);
  double num = -x.transpose() * delta;
  double denom = -x.transpose() * DavidsonSolver::dpr_correction(x, D, lambda);
  double eps = num / denom;
  delta += eps * x;
  return delta;
}

Eigen::MatrixXd DavidsonSolver::QR_ortho(const Eigen::MatrixXd &A) const {

  int nrows = A.rows();
  int ncols = A.cols();
  ncols = std::min(nrows, ncols);

  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  Eigen::MatrixXd result =
      qr.householderQ() * Eigen::MatrixXd::Identity(nrows, ncols);
  return result;
}

Eigen::MatrixXd DavidsonSolver::gramschmidt_ortho(const Eigen::MatrixXd &A,
                                                  int nstart) const {
  Eigen::MatrixXd Q = A;
  for (int j = nstart; j < A.cols(); ++j) {
    Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));
    if (Q.col(j).norm() <= 1E-12 * A.col(j).norm()) {
      throw std::runtime_error(
          "Linear dependencies in Gram-Schmidt. Switch to QR");
    }
    Q.col(j).normalize();
  }
  return Q;
}

}  // namespace xtp
}  // namespace votca