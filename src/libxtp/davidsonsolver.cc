/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <iostream>
#include <stdexcept>

// Local VOTCA includes
#include "votca/xtp/davidsonsolver.h"
#include "votca/xtp/eigen.h"

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

DavidsonSolver::DavidsonSolver(Logger &log) : log_(log) {}

void DavidsonSolver::printTiming(
    const std::chrono::time_point<std::chrono::system_clock> &start) const {
  XTP_LOG(Log::error, log_)
      << TimeStamp() << "-----------------------------------" << std::flush;
  std::chrono::time_point<std::chrono::system_clock> end =
      std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;
  XTP_LOG(Log::error, log_) << TimeStamp() << "- Davidson ran for "
                            << elapsed_time.count() << "secs." << std::flush;
  XTP_LOG(Log::error, log_)
      << TimeStamp() << "-----------------------------------" << std::flush;
}

void DavidsonSolver::checkOptions(Index operator_size) {
  //. search space exceeding the system size
  if (max_search_space_ > operator_size) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " == Warning : Max search space ("
        << max_search_space_ << ") larger than system size (" << operator_size
        << ")" << std::flush;

    max_search_space_ = operator_size;
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " == Warning : Max search space set to "
        << operator_size << std::flush;

    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " == Warning : If problems appear, try asking for less than "
        << Index(operator_size / 10) << " eigenvalues" << std::flush;
  }
}

void DavidsonSolver::printOptions(Index operator_size) const {

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Davidson Solver using " << OPENMP::getMaxThreads()
      << " threads." << std::flush;
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Tolerance : " << tol_ << std::flush;

  switch (this->davidson_correction_) {
    case CORR::DPR:
      XTP_LOG(Log::error, log_)
          << TimeStamp() << " DPR Correction" << std::flush;
      break;
    case CORR::OLSEN:
      XTP_LOG(Log::error, log_)
          << TimeStamp() << " Olsen Correction" << std::flush;
      break;
  }

  XTP_LOG(Log::error, log_) << TimeStamp() << " Matrix size : " << operator_size
                            << 'x' << operator_size << std::flush;
}

void DavidsonSolver::printIterationData(
    const DavidsonSolver::RitzEigenPair &rep,
    const DavidsonSolver::ProjectedSpace &proj, Index neigen) const {

  Index converged_roots = proj.root_converged.head(neigen).count();
  double percent_converged = 100 * double(converged_roots) / double(neigen);
  XTP_LOG(Log::error, log_)
      << TimeStamp()
      << format(" %1$4d %2$12d \t %3$4.2e \t %4$5.2f%% converged") % i_iter_ %
             proj.search_space() % rep.res_norm().head(neigen).maxCoeff() %
             percent_converged
      << std::flush;
}

void DavidsonSolver::set_matrix_type(std::string mt) {
  if (mt == "HAM") {
    this->matrix_type_ = MATRIX_TYPE::HAM;
  } else if (mt == "SYMM") {
    this->matrix_type_ = MATRIX_TYPE::SYMM;
  } else {
    throw std::runtime_error(mt + " is not a valid Davidson matrix type");
  }
}

void DavidsonSolver::set_correction(std::string method) {
  if (method == "DPR") {
    this->davidson_correction_ = CORR::DPR;
  } else if (method == "OLSEN") {
    this->davidson_correction_ = CORR::OLSEN;
  } else {
    throw std::runtime_error(method +
                             " is not a valid Davidson correction method");
  }
}

void DavidsonSolver::set_tolerance(std::string tol) {
  if (tol == "loose") {
    this->tol_ = 1E-3;
  } else if (tol == "normal") {
    this->tol_ = 1E-4;
  } else if (tol == "strict") {
    this->tol_ = 1E-5;
  } else if (tol == "lapack") {
    this->tol_ = 1E-9;
  } else {
    throw std::runtime_error(tol + " is not a valid Davidson tolerance");
  }
}

void DavidsonSolver::set_size_update(std::string update_size) {

  if (update_size == "min") {
    this->davidson_update_ = UPDATE::MIN;
  } else if (update_size == "safe") {
    this->davidson_update_ = UPDATE::SAFE;
  } else if (update_size == "max") {
    this->davidson_update_ = UPDATE::MAX;
  } else {
    throw std::runtime_error(update_size + " is not a valid Davidson update");
  }
}

Index DavidsonSolver::getSizeUpdate(Index neigen) const {
  Index size_update;
  switch (this->davidson_update_) {
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
      Eigen::MatrixXd::Zero(Adiag_.size(), size_initial_guess);
  ArrayXl idx = DavidsonSolver::argsort(Adiag_);

  switch (this->matrix_type_) {
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
      Index ind0 = Adiag_.size() / 2;
      for (Index j = 0; j < size_initial_guess; j++) {
        guess(idx(ind0 + j), j) = 1.0;
      }
      break;
  }
  return guess;
}
DavidsonSolver::RitzEigenPair DavidsonSolver::getRitzEigenPairs(
    const ProjectedSpace &proj) const {
  // get the ritz vectors
  switch (this->matrix_type_) {
    case MATRIX_TYPE::SYMM: {
      return getRitz(proj);
    }
    case MATRIX_TYPE::HAM: {
      return getHarmonicRitz(proj);
    }
  }
  return RitzEigenPair();
}

DavidsonSolver::RitzEigenPair DavidsonSolver::getRitz(
    const DavidsonSolver::ProjectedSpace &proj) const {

  DavidsonSolver::RitzEigenPair rep;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(proj.T);
  if (es.info() != Eigen::ComputationInfo::Success) {
    std::cerr << "A\n" << proj.T << std::endl;
    throw std::runtime_error("Small hermitian eigenvalue problem failed.");
  }
  // we only need enough pairs for either extension of space or restart
  Index needed_pairs =
      std::min(proj.T.cols(), std::max(restart_size_, proj.size_update));
  rep.lambda = es.eigenvalues().head(needed_pairs);
  rep.U = es.eigenvectors().leftCols(needed_pairs);

  rep.q = proj.V * rep.U;                                       // Ritz vectors
  rep.res = proj.AV * rep.U - rep.q * rep.lambda.asDiagonal();  // residues
  return rep;
}

DavidsonSolver::RitzEigenPair DavidsonSolver::getHarmonicRitz(
    const ProjectedSpace &proj) const {

  /* Compute the Harmonic Ritz vector following
   * Computing Interior Eigenvalues of Large Matrices
   * Ronald B Morgan
   * LINEAR ALGEBRA AND ITS APPLICATIONS 154-156:289-309 (1991)
   * https://cpb-us-w2.wpmucdn.com/sites.baylor.edu/dist/e/71/files/2015/05/InterEvals-1vgdz91.pdf
   */

  RitzEigenPair rep;
  bool return_eigenvectors = true;
  Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges(proj.T, proj.B,
                                                     return_eigenvectors);
  if (ges.info() != Eigen::ComputationInfo::Success) {
    std::cerr << "A\n" << proj.T << std::endl;
    std::cerr << "B\n" << proj.B << std::endl;
    throw std::runtime_error("Small generalized eigenvalue problem failed.");
  }

  std::vector<std::pair<Index, Index>> complex_pairs;
  for (Index i = 0; i < ges.eigenvalues().size(); i++) {
    if (ges.eigenvalues()(i).imag() != 0) {
      bool found_partner = false;
      for (auto &pair : complex_pairs) {
        if (pair.second > -1) {
          continue;
        } else {
          bool are_pair = (std::abs(ges.eigenvalues()(pair.first).real() -
                                    ges.eigenvalues()(i).real()) < 1e-9) &&
                          (std::abs(ges.eigenvalues()(pair.first).imag() +
                                    ges.eigenvalues()(i).imag()) < 1e-9);
          if (are_pair) {
            pair.second = i;
            found_partner = true;
          }
        }
      }

      if (!found_partner) {
        complex_pairs.emplace_back(i, -1);
      }
    }
  }

  for (const auto &pair : complex_pairs) {
    if (pair.second < 0) {
      throw std::runtime_error("Eigenvalue:" + std::to_string(pair.first) +
                               " is complex but has no partner.");
    }
  }
  if (!complex_pairs.empty()) {
    XTP_LOG(Log::warning, log_)
        << TimeStamp() << " Found " << complex_pairs.size()
        << " complex pairs in eigenvalue problem" << std::flush;
  }
  Eigen::VectorXd eigenvalues =
      Eigen::VectorXd::Zero(ges.eigenvalues().size() - complex_pairs.size());
  Eigen::MatrixXd eigenvectors =
      Eigen::MatrixXd::Zero(ges.eigenvectors().rows(),
                            ges.eigenvectors().cols() - complex_pairs.size());

  Index j = 0;
  for (Index i = 0; i < ges.eigenvalues().size(); i++) {
    bool is_second_in_complex_pair =
        std::find_if(complex_pairs.begin(), complex_pairs.end(),
                     [&](const std::pair<Index, Index> &pair) {
                       return pair.second == i;
                     }) != complex_pairs.end();
    if (is_second_in_complex_pair) {
      continue;
    } else {
      eigenvalues(j) = ges.eigenvalues()(i).real();
      eigenvectors.col(j) = ges.eigenvectors().col(i).real();
      eigenvectors.col(j).normalize();
      j++;
    }
  }
  // we only need enough pairs for either extension of space or restart
  Index needed_pairs =
      std::min(proj.T.cols(), std::max(restart_size_, proj.size_update));
  ArrayXl idx =
      DavidsonSolver::argsort(eigenvalues).reverse().head(needed_pairs);
  // we need the largest values, because this is the inverse value, so
  // reverse list

  rep.U = DavidsonSolver::extract_vectors(eigenvectors, idx);
  rep.lambda = (rep.U.transpose() * proj.T * rep.U).diagonal();
  rep.q = proj.V * rep.U;                                       // Ritz vectors
  rep.res = proj.AV * rep.U - rep.q * rep.lambda.asDiagonal();  // residues
  return rep;
}

DavidsonSolver::ProjectedSpace DavidsonSolver::initProjectedSpace(
    Index neigen, Index size_initial_guess) const {
  DavidsonSolver::ProjectedSpace proj;

  // initial vector basis
  proj.V = DavidsonSolver::setupInitialEigenvectors(size_initial_guess);

  // update variables
  proj.size_update = DavidsonSolver::getSizeUpdate(neigen);
  proj.root_converged = ArrayXb::Constant(proj.size_update, false);
  return proj;
}

bool DavidsonSolver::checkConvergence(const DavidsonSolver::RitzEigenPair &rep,
                                      DavidsonSolver::ProjectedSpace &proj,
                                      Index neigen) const {
  proj.root_converged = (rep.res_norm().head(proj.size_update) < tol_);
  return proj.root_converged.head(neigen).all();
}

Index DavidsonSolver::extendProjection(
    const DavidsonSolver::RitzEigenPair &rep,
    DavidsonSolver::ProjectedSpace &proj) const {

  Index nupdate = (proj.root_converged == false).count();
  Index oldsize = proj.V.cols();
  proj.V.conservativeResize(Eigen::NoChange, oldsize + nupdate);

  Index k = 0;
  for (Index j = 0; j < proj.size_update; j++) {
    // skip the roots that have already converged
    if (proj.root_converged[j]) {
      continue;
    }
    // residue vector
    Eigen::VectorXd w =
        computeCorrectionVector(rep.q.col(j), rep.lambda(j), rep.res.col(j));

    // append the correction vector to the search space
    proj.V.col(oldsize + k) = w.normalized();
    k++;
  }
  orthogonalize(proj.V, nupdate);
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

  /* compute correction vector with either DPR or OLSEN CORRECTION
   * For details on the method see :
   * Systematic Study of Selected Diagonalization Methods
   * for Configuration Interaction Matrices
   * M.L. Leininger et al .
   * Journal of Computational Chemistry Vol 22, No. 13 1574-1589 (2001)
   */
  Eigen::VectorXd correction;
  switch (this->davidson_correction_) {
    case CORR::DPR: {
      correction = dpr(Aqj, lambdaj);
      break;
    }
    case CORR::OLSEN: {
      correction = olsen(Aqj, qj, lambdaj);
      break;
    }
  }
  // make sure no nan values are there, instead we set them to zero
  return correction.unaryExpr(
      [](double v) { return std::isfinite(v) ? v : 0.0; });
}

Eigen::VectorXd DavidsonSolver::dpr(const Eigen::VectorXd &r,
                                    double lambda) const {
  /* \brief Compute the diagonal preconditoned residue :
  \delta = -r/(D - lambda)
   */
  return (-r.array() / (Adiag_.array() - lambda));
}

Eigen::VectorXd DavidsonSolver::olsen(const Eigen::VectorXd &r,
                                      const Eigen::VectorXd &x,
                                      double lambda) const {
  /* \brief Compute the olsen correction :
  \delta = (D-\lambda)^{-1} (-r + \epsilon x)
  */
  Eigen::VectorXd delta = DavidsonSolver::dpr(r, lambda);
  double num = -x.transpose() * delta;
  double denom = -x.transpose() * dpr(x, lambda);
  double eps = num / denom;
  delta += eps * x;
  return delta;
}

void DavidsonSolver::orthogonalize(Eigen::MatrixXd &V, Index nupdate) const {
  DavidsonSolver::gramschmidt(V, V.cols() - nupdate);
}

void DavidsonSolver::gramschmidt(Eigen::MatrixXd &Q, Index nstart) const {
  Index nupdate = Q.cols() - nstart;
  Eigen::VectorXd norms = Q.rightCols(nupdate).colwise().norm();
  // orthogonalize with respect to already existing vectors
  if (nstart > 0) {
    Q.rightCols(nupdate) -=
        Q.leftCols(nstart) *
        (Q.leftCols(nstart).transpose() * Q.rightCols(nupdate));
    Q.rightCols(nupdate).colwise().normalize();
  }
  // orthogonalize vectors to each other
  for (Index j = nstart + 1; j < Q.cols(); ++j) {
    Index range = j - nstart;
    Q.col(j) -= Q.block(0, nstart, Q.rows(), range) *
                (Q.block(0, nstart, Q.rows(), range).transpose() * Q.col(j));
    Q.col(j).normalize();
  }
  // repeat again two is enough GS
  // http://stoppels.blog/posts/orthogonalization-performance
  if (nstart > 0) {
    Q.rightCols(nupdate) -=
        Q.leftCols(nstart) *
        (Q.leftCols(nstart).transpose() * Q.rightCols(nupdate));
    Q.rightCols(nupdate).colwise().normalize();
  }

  for (Index j = nstart + 1; j < Q.cols(); ++j) {
    Index range = j - nstart;
    Q.col(j) -= Q.block(0, nstart, Q.rows(), range) *
                (Q.block(0, nstart, Q.rows(), range).transpose() * Q.col(j));
    if (Q.col(j).norm() <= 1E-12 * norms(range)) {
      // info_ = Eigen::ComputationInfo::NumericalIssue;
      throw std::runtime_error("Linear dependencies in Gram-Schmidt.");
    }
    Q.col(j).normalize();
  }
}

Eigen::MatrixXd DavidsonSolver::qr(const Eigen::MatrixXd &A) const {

  Index nrows = A.rows();
  Index ncols = A.cols();
  ncols = std::min(nrows, ncols);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nrows, ncols);
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  return qr.householderQ() * I;
}

void DavidsonSolver::restart(const DavidsonSolver::RitzEigenPair &rep,
                             DavidsonSolver::ProjectedSpace &proj,
                             Index newvectors) const {
  Eigen::MatrixXd newV =
      Eigen::MatrixXd(proj.V.rows(), newvectors + restart_size_);
  newV.rightCols(newvectors) = proj.V.rightCols(newvectors);
  if (matrix_type_ == MATRIX_TYPE::SYMM) {

    newV.leftCols(restart_size_) = rep.q.leftCols(restart_size_);
    proj.AV *= rep.U.leftCols(restart_size_);  // corresponds to replacing
                                               // V with q.leftCols
  } else {
    Eigen::MatrixXd orthonormal =
        DavidsonSolver::qr(rep.U.leftCols(restart_size_));
    newV.leftCols(restart_size_) =
        proj.V.leftCols(proj.V.cols() - newvectors) * orthonormal;
    proj.AV *= orthonormal;

    proj.AAV *= orthonormal;
    proj.B = newV.leftCols(restart_size_).transpose() * proj.AAV;
  }
  proj.T = newV.leftCols(restart_size_).transpose() * proj.AV;
  proj.V = newV;
}

void DavidsonSolver::storeConvergedData(
    const DavidsonSolver::RitzEigenPair &rep, Index neigen) {

  DavidsonSolver::storeEigenPairs(rep, neigen);
  XTP_LOG(Log::error, log_) << TimeStamp() << " Davidson converged after "
                            << i_iter_ << " iterations." << std::flush;
  info_ = Eigen::ComputationInfo::Success;
}

void DavidsonSolver::storeNotConvergedData(
    const DavidsonSolver::RitzEigenPair &rep, const ArrayXb &root_converged,
    Index neigen) {

  DavidsonSolver::storeEigenPairs(rep, neigen);

  double percent_converged = 0;

  for (Index i = 0; i < neigen; i++) {
    if (!root_converged[i]) {
      eigenvalues_(i) = 0;
      eigenvectors_.col(i).setZero();
    } else {
      percent_converged += 1.;
    }
  }
  percent_converged /= double(neigen);
  percent_converged *= 100.;
  XTP_LOG(Log::error, log_)
      << TimeStamp() << "- Warning : Davidson "
      << format("%1$5.2f%%") % percent_converged << " converged after "
      << i_iter_ << " iterations." << std::flush;
  info_ = Eigen::ComputationInfo::NoConvergence;
}

void DavidsonSolver::storeEigenPairs(const DavidsonSolver::RitzEigenPair &rep,
                                     Index neigen) {
  // store the eigenvalues/eigenvectors
  this->eigenvalues_ = rep.lambda.head(neigen);
  this->eigenvectors_ = rep.q.leftCols(neigen);
  this->eigenvectors_.colwise().normalize();
}

}  // namespace xtp
}  // namespace votca
