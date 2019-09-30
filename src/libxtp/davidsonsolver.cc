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
#include <Spectra/SymEigsSolver.h>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

using namespace std;

DavidsonSolver::DavidsonSolver(Logger &log) : _log(log) {}

Eigen::VectorXd DavidsonSolver::ComputeCorrectionVector(
    const Eigen::VectorXd &Adiag, const Eigen::VectorXd &qj, double lambdaj,
    const Eigen::VectorXd &Aqj) const {
  Eigen::VectorXd w;
  // compute correction vector
  switch (this->_davidson_correction) {
    case CORR::DPR:
      w = dpr_correction(Aqj, Adiag, lambdaj);
      break;
    case CORR::OLSEN:
      w = olsen_correction(Aqj, qj, Adiag, lambdaj);
      break;
  }
  return w;
}

void DavidsonSolver::PrintTiming(
    const std::chrono::time_point<std::chrono::system_clock> &start) const {
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << "-----------------------------------" << flush;
  std::chrono::time_point<std::chrono::system_clock> end =
      std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;
  XTP_LOG_SAVE(logDEBUG, _log) << TimeStamp() << "- Davidson ran for "
                               << elapsed_time.count() << "secs." << flush;
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << "-----------------------------------" << flush;
}

void DavidsonSolver::PrintOptions(int op_size) const {

  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << " Davidson Solver using " << OPENMP::getMaxThreads()
      << " threads." << flush;
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << " Tolerance : " << _tol << flush;

  switch (this->_davidson_correction) {
    case CORR::DPR:
      XTP_LOG_SAVE(logDEBUG, _log) << TimeStamp() << " DPR Correction" << flush;
      break;
    case CORR::OLSEN:
      XTP_LOG_SAVE(logDEBUG, _log)
          << TimeStamp() << " Olsen Correction" << flush;
      break;
  }

  switch (this->_davidson_ortho) {
    case ORTHO::GS:
      XTP_LOG_SAVE(logDEBUG, _log)
          << TimeStamp() << " Gram-Schmidt Orthogonalization" << flush;
      break;
    case ORTHO::QR:
      XTP_LOG_SAVE(logDEBUG, _log)
          << TimeStamp() << " QR Orthogonalization" << flush;
      break;
  }
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << " Matrix size : " << op_size << 'x' << op_size << flush;
}

void DavidsonSolver::set_matrix_type(std::string mt) {
  if (mt == "HAM")
    this->_matrix_type = TYPE::HAM;
  else if (mt == "SYMM")
    this->_matrix_type = TYPE::SYMM;
  else
    throw std::runtime_error(
        mt + " is not a valid Davidson matrix type");
}

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

int DavidsonSolver::get_size_update(int neigen) const {
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
    default:
      size_update = 2 * neigen;
      break;
  }
  return size_update;
}

Eigen::ArrayXi DavidsonSolver::argsort(const Eigen::VectorXd &V) const {
  /* \brief return the index of the sorted vector */
  Eigen::ArrayXi idx = Eigen::ArrayXi::LinSpaced(V.rows(), 0, V.rows() - 1);
  std::sort(idx.data(), idx.data() + idx.size(),
            [&](int i1, int i2) { return V[i1] < V[i2]; });
  return idx;
}

Eigen::MatrixXd DavidsonSolver::SetupInitialEigenvectors(
    Eigen::VectorXd &d, int size_initial_guess) const {

  Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(d.size(), size_initial_guess);
  Eigen::ArrayXi idx = DavidsonSolver::argsort(d);

  switch (this->_matrix_type) {
    case TYPE::SYMM:
      /* \brief Initialize the guess eigenvector so that they 'target' the lowest
       * diagonal elements */
      for (int j = 0; j < size_initial_guess; j++) {
        guess(idx(j), j) = 1.0;
      }
      break;

    case TYPE::HAM:
      /* Initialize the guess eigenvector so that they 'target' the lowest
       * positive diagonal elements */
        int ind0 = d.size()/2;
        int shift = size_initial_guess/4;
      for (int j = 0; j < size_initial_guess; j++) {
        guess(idx(ind0-shift+j), j) = 1.0;
      }
      break;
    }
  return guess;
}

Eigen::ArrayXi DavidsonSolver::index_window(const Eigen::VectorXd &V, 
    int size_update, double target_min_val, double perc_below) const {

  double min_val = 1E12;
  int index_min = -1;
  int n = V.size();
  Eigen::ArrayXi idx = Eigen::ArrayXi::Zero(size_update);

  // reorder values in ascending order
  Eigen::ArrayXi isort = argsort(V);

  // index of lowest element closed to 0
  for (int i=0; i < n; i++) {
    if ( (V(isort(i)) > target_min_val) && ( V(isort(i)) < min_val) ) {
      min_val = V(isort(i));
      index_min = i;
    }
  }

  // index of the element we want in size update
  int shift = size_update * perc_below;
  int index_start = index_min - shift;

  if (index_start < 0)
    index_start = 0;

  for (int i=0; i<size_update; i++) {
    idx(i) = isort(index_start+i);
  }

  return idx;
}

Eigen::MatrixXd DavidsonSolver::extract_eigenvectors(const Eigen::MatrixXd &V, 
    const Eigen::ArrayXi &idx) const {
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(V.rows(),idx.size());
  for (int i=0; i < idx.size(); i++) {
    W.col(i) = V.col(idx(i));
  }
  return W;
}

Eigen::VectorXd DavidsonSolver::dpr_correction(const Eigen::VectorXd &r,
                                               const Eigen::VectorXd &D,
                                               double lambda) const {
  /* \brief Compute the diagonal preconditoned residue : delta = - (D -
   * lambda)^{-1} r
   */
  Eigen::VectorXd delta = r.array() / (lambda - D.array());
  return delta;
}

Eigen::VectorXd DavidsonSolver::olsen_correction(const Eigen::VectorXd &r,
                                                 const Eigen::VectorXd &x,
                                                 const Eigen::VectorXd &D,
                                                 double lambda) const {
  /* \brief Compute the olsen correction :

  \delta = (D-\lambda)^{-1} (-r + \epsilon x)

  */
  int size = r.rows();
  Eigen::VectorXd delta = Eigen::VectorXd::Zero(size);
  delta = DavidsonSolver::dpr_correction(r, D, lambda);
  double num = -x.transpose() * delta;
  double denom = -x.transpose() * dpr_correction(x, D, lambda);
  double eps = num / denom;
  delta += eps * x;
  return delta;
}

Eigen::MatrixXd DavidsonSolver::QR_ortho(const Eigen::MatrixXd &A) const {

  int nrows = A.rows();
  int ncols = A.cols();
  ncols = std::min(nrows, ncols);

  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  Eigen::MatrixXd result = qr.householderQ();
  result.conservativeResize(nrows, ncols);
  return result;
}

Eigen::MatrixXd DavidsonSolver::gramschmidt_ortho(const Eigen::MatrixXd &A,
                                                  int nstart) {
  Eigen::MatrixXd Q = A;
  for (int j = nstart; j < A.cols(); ++j) {
    Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));
    if (Q.col(j).norm() <= 1E-12 * A.col(j).norm()) {
      _info = Eigen::ComputationInfo::NumericalIssue;
      throw std::runtime_error(
          "Linear dependencies in Gram-Schmidt. Switch to QR");
    }
    Q.col(j).normalize();
  }
  return Q;
}

}  // namespace xtp
}  // namespace votca