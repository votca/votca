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

DavidsonSolver::DavidsonSolver(Logger &log) : _log(log) {}

void DavidsonSolver::printTiming(
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

void DavidsonSolver::checkOptions(Index operator_size) {
  //. search space exceeding the system size
  if (_max_search_space > operator_size) {
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " == Warning : Max search space ("
        << _max_search_space << ") larger than system size (" << operator_size
        << ")" << flush;

    _max_search_space = operator_size;
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " == Warning : Max search space set to "
        << operator_size << flush;

    this->_davidson_ortho = ORTHO::QR;
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp()
        << " == Warning : Orthogonalization set to QR for stabilty " << flush;

    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp()
        << " == Warning : If problems appear, try asking for less than "
        << Index(operator_size / 10) << " eigenvalues" << flush;
  }
}

void DavidsonSolver::printOptions(Index operator_size) const {

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
      << TimeStamp() << " Matrix size : " << operator_size << 'x'
      << operator_size << flush;
}

void DavidsonSolver::printIterationData(
    const DavidsonSolver::RitzEigenPair &rep,
    const DavidsonSolver::ProjectedSpace &proj, Index neigen,
    Index iiter) const {

  if (iiter == 0) {
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " iter\tSearch Space\tNorm" << flush;
  }

  Index converged_roots = 0;
  for (Index i = 0; i < neigen; i++) {
    converged_roots += proj.root_converged[i];
  }
  double percent_converged = 100 * double(converged_roots) / double(neigen);
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp()
      << format(" %1$4d %2$12d \t %3$4.2e \t %4$5.2f%% converged") % iiter %
             proj.search_space % rep.res_norm.head(neigen).maxCoeff() %
             percent_converged
      << flush;
}

void DavidsonSolver::set_matrix_type(std::string mt) {
  if (mt == "HAM") {
    this->_matrix_type = MATRIX_TYPE::HAM;
  } else if (mt == "SYMM") {
    this->_matrix_type = MATRIX_TYPE::SYMM;
  } else {
    throw std::runtime_error(mt + " is not a valid Davidson matrix type");
  }
}

void DavidsonSolver::set_ortho(std::string method) {
  if (method == "GS") {
    this->_davidson_ortho = ORTHO::GS;
  } else if (method == "QR") {
    this->_davidson_ortho = ORTHO::QR;
  } else {
    throw std::runtime_error(
        method + " is not a valid Davidson orthogonalization method");
  }
}

void DavidsonSolver::set_correction(std::string method) {
  if (method == "DPR") {
    this->_davidson_correction = CORR::DPR;
  } else if (method == "OLSEN") {
    this->_davidson_correction = CORR::OLSEN;
  } else {
    throw std::runtime_error(method +
                             " is not a valid Davidson correction method");
  }
}

void DavidsonSolver::set_tolerance(std::string tol) {
  if (tol == "loose") {
    this->_tol = 1E-3;
  } else if (tol == "normal") {
    this->_tol = 1E-4;
  } else if (tol == "strict") {
    this->_tol = 1E-5;
  } else {
    throw std::runtime_error(tol + " is not a valid Davidson tolerance");
  }
}

void DavidsonSolver::set_size_update(std::string update_size) {

  if (update_size == "min") {
    this->_davidson_update = UPDATE::MIN;
  } else if (update_size == "safe") {
    this->_davidson_update = UPDATE::SAFE;
  } else if (update_size == "max") {
    this->_davidson_update = UPDATE::MAX;
  } else {
    throw std::runtime_error(update_size + " is not a valid Davidson update");
  }
}

long DavidsonSolver::getSizeUpdate(Index neigen) const {
  Index size_update;
  switch (this->_davidson_update) {
    case UPDATE::MIN:
      size_update = neigen;
      break;
    case UPDATE::SAFE:
      if (neigen < 20) {
        size_update = static_cast<Index>(1.5 * double(neigen));
      } else {
        size_update = neigen + 10;
      }
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

ArrayXl DavidsonSolver::argsort(const Eigen::VectorXd &V) const {
  /* \brief return the index of the sorted vector */
  ArrayXl idx = ArrayXl::LinSpaced(V.rows(), 0, V.rows() - 1);
  std::sort(idx.data(), idx.data() + idx.size(),
            [&](Index i1, Index i2) { return V[i1] < V[i2]; });
  return idx;
}

Eigen::MatrixXd DavidsonSolver::setupInitialEigenvectors(
    Index size_initial_guess) const {

  Eigen::MatrixXd guess =
      Eigen::MatrixXd::Zero(_Adiag.size(), size_initial_guess);
  ArrayXl idx = DavidsonSolver::argsort(_Adiag);

  switch (this->_matrix_type) {
    case MATRIX_TYPE::SYMM:
      /* \brief Initialize the guess eigenvector so that they 'target' the
       * smallest diagonal elements */
      for (Index j = 0; j < size_initial_guess; j++) {
        guess(idx(j), j) = 1.0;
      }
      break;

    case MATRIX_TYPE::HAM:
      /* Initialize the guess eigenvector so that they 'target' the lowest
       * positive diagonal elements */
      Index ind0 = _Adiag.size() / 2;
      for (Index j = 0; j < size_initial_guess; j++) {
        guess(idx(ind0 + j), j) = 1.0;
      }
      break;
  }
  return guess;
}

DavidsonSolver::RitzEigenPair DavidsonSolver::getRitz(
    const DavidsonSolver::ProjectedSpace &proj) const {

  DavidsonSolver::RitzEigenPair rep;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(proj.T);
  rep.lambda = es.eigenvalues();
  rep.U = es.eigenvectors();

  rep.q = proj.V * rep.U;                                       // Ritz vectors
  rep.res = proj.AV * rep.U - rep.q * rep.lambda.asDiagonal();  // residues
  rep.res_norm = rep.res.colwise().norm();  // residues norms

  return rep;
}

DavidsonSolver::ProjectedSpace DavidsonSolver::initProjectedSpace(
    Index neigen, Index size_initial_guess) const {
  DavidsonSolver::ProjectedSpace proj;

  // initial vector basis
  proj.V = DavidsonSolver::setupInitialEigenvectors(size_initial_guess);
  proj.search_space = proj.V.cols();

  // update variables
  proj.size_update = DavidsonSolver::getSizeUpdate(neigen);
  proj.root_converged = std::vector<bool>(proj.size_update, false);

  return proj;
}

long DavidsonSolver::extendProjection(DavidsonSolver::RitzEigenPair &rep,
                                      DavidsonSolver::ProjectedSpace &proj) {

  Index nupdate = 0;
  for (Index j = 0; j < proj.size_update; j++) {

    // skip the root that have already converged
    if (this->_matrix_type == MATRIX_TYPE::SYMM) {
      if (proj.root_converged[j]) {
        continue;
      }
    }
    nupdate++;

    // residue vector
    Eigen::VectorXd w =
        computeCorrectionVector(rep.q.col(j), rep.lambda(j), rep.res.col(j));

    // append the correction vector to the search space
    proj.V.conservativeResize(Eigen::NoChange, proj.V.cols() + 1);
    proj.V.rightCols<1>() = w.normalized();

    // track converged root
    proj.root_converged[j] = (rep.res_norm[j] < _tol);
  }

  proj.search_space = proj.V.cols();

  return nupdate;
}

Eigen::MatrixXd DavidsonSolver::extract_vectors(const Eigen::MatrixXd &V,
                                                const ArrayXl &idx) const {
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(V.rows(), idx.size());
  for (Index i = 0; i < idx.size(); i++) {
    W.col(i) = V.col(idx(i));
  }
  return W;
}

Eigen::VectorXd DavidsonSolver::computeCorrectionVector(
    const Eigen::VectorXd &qj, double lambdaj,
    const Eigen::VectorXd &Aqj) const {

  /* compute correction vector with either DPR or OLSEn CORRECTION
   * For details on the method see :
   * Systematic Study of Selected Diagonalization Methods
   * for Configuration Interaction Matrices
   * M.L. Leininger et al .
   * Journal of Computational Chemistry Vol 22, No. 13 1574-1589 (2001)
   */

  Eigen::VectorXd out_vect;
  switch (this->_davidson_correction) {
    case CORR::DPR: {
      out_vect = dpr(Aqj, lambdaj);
      break;
    }

    case CORR::OLSEN: {
      out_vect = olsen(Aqj, qj, lambdaj);
      break;
    }
  }
  return out_vect;
}

Eigen::VectorXd DavidsonSolver::dpr(const Eigen::VectorXd &r,
                                    double lambda) const {
  /* \brief Compute the diagonal preconditoned residue : delta = - (D -
   * lambda)^{-1} r
   */
  Eigen::VectorXd delta = r.array() / (lambda - _Adiag.array());
  return delta;
}

Eigen::VectorXd DavidsonSolver::olsen(const Eigen::VectorXd &r,
                                      const Eigen::VectorXd &x,
                                      double lambda) const {
  /* \brief Compute the olsen correction :

  \delta = (D-\lambda)^{-1} (-r + \epsilon x)

  */
  Index size = r.rows();
  Eigen::VectorXd delta = Eigen::VectorXd::Zero(size);
  delta = DavidsonSolver::dpr(r, lambda);
  double num = -x.transpose() * delta;
  double denom = -x.transpose() * dpr(x, lambda);
  double eps = num / denom;
  delta += eps * x;
  return delta;
}

Eigen::MatrixXd DavidsonSolver::orthogonalize(const Eigen::MatrixXd &V,
                                              Index nupdate) {

  Eigen::MatrixXd Vout;
  switch (this->_davidson_ortho) {
    case ORTHO::GS: {
      Vout = DavidsonSolver::gramschmidt(V, V.cols() - nupdate);
      break;
    }
    case ORTHO::QR: {
      Vout = DavidsonSolver::qr(V);
      break;
    }
  }
  return Vout;
}

Eigen::MatrixXd DavidsonSolver::qr(const Eigen::MatrixXd &A) const {

  Index nrows = A.rows();
  Index ncols = A.cols();
  ncols = std::min(nrows, ncols);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nrows, ncols);
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  return qr.householderQ() * I;
}

Eigen::MatrixXd DavidsonSolver::gramschmidt(const Eigen::MatrixXd &A,
                                            Index nstart) {
  Eigen::MatrixXd Q = A;
  for (Index j = nstart; j < A.cols(); ++j) {
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

void DavidsonSolver::restart(const DavidsonSolver::RitzEigenPair &rep,
                             DavidsonSolver::ProjectedSpace &proj,
                             Index size_restart) const {
  proj.V = rep.q.leftCols(size_restart);
  proj.V.colwise().normalize();
  proj.AV = proj.AV * rep.U.leftCols(size_restart);  // corresponds to replacing
                                                     // V with q.leftCols
  proj.T = proj.V.transpose() * proj.AV;
  proj.search_space = size_restart;
}

void DavidsonSolver::storeConvergedData(
    const DavidsonSolver::RitzEigenPair &rep, Index neigen, Index iiter) {

  DavidsonSolver::storeEigenPairs(rep, neigen);
  this->_num_iter = iiter;
  XTP_LOG_SAVE(logDEBUG, _log) << TimeStamp() << " Davidson converged after "
                               << iiter << " iterations." << flush;
  _info = Eigen::ComputationInfo::Success;
}

void DavidsonSolver::storeNotConvergedData(
    const DavidsonSolver::RitzEigenPair &rep, std::vector<bool> &root_converged,
    Index neigen) {

  DavidsonSolver::storeEigenPairs(rep, neigen);
  this->_num_iter = _iter_max;

  double percent_converged = 0;

  for (Index i = 0; i < neigen; i++) {
    if (!root_converged[i]) {
      _eigenvalues(i) = 0;
      _eigenvectors.col(i).setZero();
    } else {
      percent_converged += 1.;
    }
  }
  percent_converged /= double(neigen);
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << "- Warning : Davidson "
      << format("%1$5.2f%%") % percent_converged << " converged after "
      << _iter_max << " iterations." << flush;
  _info = Eigen::ComputationInfo::NoConvergence;
}

void DavidsonSolver::storeEigenPairs(const DavidsonSolver::RitzEigenPair &rep,
                                     Index neigen) {
  // store the eigenvalues/eigenvectors
  this->_eigenvalues = rep.lambda.head(neigen);
  this->_eigenvectors = rep.q.leftCols(neigen);
  this->_eigenvectors.colwise().normalize();
  this->_res = rep.res_norm.head(neigen);
}

}  // namespace xtp
}  // namespace votca