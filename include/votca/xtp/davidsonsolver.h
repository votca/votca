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

#ifndef __VOTCA_TOOLS_DAVIDSON_SOLVER_H
#define __VOTCA_TOOLS_DAVIDSON_SOLVER_H

#include <chrono>
#include <iostream>
#include <stdexcept>

#include <votca/xtp/eigen.h>

#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/threecenter.h>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

/**
* \brief Use Davidson algorithm to solve A*V=E*V

**/

class DavidsonSolver {

 public:
  DavidsonSolver(ctp::Logger &log);

  void set_iter_max(int N) { this->iter_max = N; }

  void set_max_search_space(int N) { this->max_search_space = N; }

  void set_tolerance(std::string tol);
  void set_correction(std::string method);
  void set_ortho(std::string method);
  void set_size_update(std::string method);
  int get_size_update(int neigen);

  Eigen::VectorXd eigenvalues() const { return this->_eigenvalues; }
  Eigen::MatrixXd eigenvectors() const { return this->_eigenvectors; }

  template <typename MatrixReplacement>
  void solve(MatrixReplacement &A, int neigen, int size_initial_guess = 0) {

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Davidson Solver" << flush;

    switch (this->davidson_correction) {

      case CORR::DPR:
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " DPR Correction" << flush;
        break;

      case CORR::OLSEN:
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " Olsen Correction" << flush;
        break;
    }

    switch (this->davidson_ortho) {

      case ORTHO::GS:
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " Gram-Schmidt Ortgonalization" << flush;
        break;

      case ORTHO::QR:
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " QR Ortgonalization" << flush;
        break;
    }

    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Tolerance : " << this->tol << flush;

    int size = A.rows();

    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Matrix size : " << size << 'x' << size
        << flush;

    //. search space exeeding the system size
    if (max_search_space > size) {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << " Warning Max search space ("
          << max_search_space << ") larger than system size (" << size << ")"
          << flush;
      max_search_space = size;
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << " Max search space set to " << size << flush;
    }

    // initial guess size
    if (size_initial_guess == 0) {
      size_initial_guess = 2 * neigen;
    }
    int search_space = size_initial_guess;
    int size_restart = size_initial_guess;
    int size_update = DavidsonSolver::get_size_update(neigen);
    int nupdate;

    Eigen::ArrayXd res_norm = Eigen::ArrayXd::Zero(size_update);
    Eigen::ArrayXd root_converged = Eigen::ArrayXd::Zero(size_update);

    double percent_converged;
    double conv;
    bool has_converged = false;

    // initialize the guess eigenvector
    Eigen::VectorXd Adiag = A.diagonal();

    // target the lowest diagonal element
    Eigen::MatrixXd V =
        DavidsonSolver::get_initial_eigenvectors(Adiag, size_initial_guess);

    // eigenvalues holder
    Eigen::VectorXd lambda;

    // temp varialbes
    Eigen::MatrixXd T, U, q;
    Eigen::VectorXd w, tmp;

    // project the matrix on the trial subspace
    T = V.transpose() * (A * V);

    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " iter\tSearch Space\tNorm" << flush;

    std::chrono::time_point<std::chrono::system_clock> istart, iend;
    std::chrono::duration<double> elapsed_time;

    for (int iiter = 0; iiter < iter_max; iiter++) {

      istart = std::chrono::system_clock::now();

      // diagonalize the small subspace
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
      lambda = es.eigenvalues();
      U = es.eigenvectors();


      // Ritz eigenvectors
      q = V.block(0, 0, V.rows(), search_space) * U;

      // correction vectors
      nupdate = 0;
      for (int j = 0; j < size_update; j++) {

        if (root_converged[j]) continue;
        nupdate += 1;

        // residue vector
        w = A * q.col(j) - lambda(j) * q.col(j);
        res_norm[j] = w.norm();

        switch (this->davidson_correction) {

          case CORR::DPR:
            w = DavidsonSolver::dpr_correction(w, Adiag, lambda(j));
            break;

          case CORR::OLSEN:
            tmp = q.col(j);
            w = DavidsonSolver::olsen_correction(w, tmp, Adiag, lambda(j));
            break;
        }

        // append the correction vector to the search space
        V.conservativeResize(Eigen::NoChange, V.cols() + 1);
        V.col(V.cols() - 1) = w.normalized();

        // track converged root
        root_converged[j] = res_norm[j] < tol;
      }

      iend = std::chrono::system_clock::now();
      elapsed_time = iend - istart;

      // Print iteration data
      percent_converged = 100 * root_converged.head(neigen).sum() / neigen;
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp()
          << format(
                 " %1$4d %2$12d \t %3$4.2e \t %4$5.2f%% converged "
                 "%5$f sec.") %
                 iiter % search_space % res_norm.head(neigen).maxCoeff() %
                 percent_converged % elapsed_time.count()
          << flush;

      // update
      search_space = V.cols();

      // break if converged
      if ((res_norm.head(neigen) < tol).all()) {
        has_converged = true;
        break;
      }

      // check if we need to restart
      if (search_space > max_search_space or search_space > size) {

        V = q.leftCols(size_restart);
        for (int j = 0; j < size_restart; j++) {
          V.col(j).normalize();
        }
        search_space = size_restart;

        // recompute the projected matrix
        T = V.transpose() * (A * V);
      }

      // continue otherwise
      else {

        switch (this->davidson_ortho){
          case ORTHO::GS:
            // orthogonalize the V vectors
            V = DavidsonSolver::gramschmidt_ortho(V, V.cols() - nupdate);
            DavidsonSolver::update_projected_matrix<MatrixReplacement>(T, A, V);
            break;
          case ORTHO::QR:
            V = DavidsonSolver::QR_ortho(V);
            T = V.transpose() * (A * V);
            break;
        }
      }
    }

    end = std::chrono::system_clock::now();
    elapsed_time = end - start;

    // store the eigenvalues/eigenvectors
    this->_eigenvalues = lambda.head(neigen);
    this->_eigenvectors = q.leftCols(neigen);

    // normalize the eigenvectors
    for (int i = 0; i < neigen; i++) {
      this->_eigenvectors.col(i).normalize();
    }

    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << "-----------------------------------" << flush;

    if (!has_converged) {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << "- Warning : Davidson " << percent_converged
          << "% converged after " << iter_max << " iterations." << flush;

      for (int i = 0; i < neigen; i++) {
        if (not root_converged[i]) {
          this->_eigenvalues(i) = 0;
          this->_eigenvectors.col(i) = Eigen::VectorXd::Zero(size);
        }
      }
    } else {

      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << "- Davidson converged in "
          << elapsed_time.count() << "secs." << flush;
    }
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << "-----------------------------------" << flush;
  }

 private:
  ctp::Logger &_log;
  int iter_max = 50;
  double tol = 1E-4;
  int max_search_space = 1000;

  enum CORR { DPR, OLSEN };
  CORR davidson_correction = CORR::DPR;

  enum UPDATE { MIN, SAFE, MAX };
  UPDATE davidson_update = UPDATE::SAFE;

  enum ORTHO { GS, QR};
  ORTHO davidson_ortho = ORTHO::GS;

  Eigen::VectorXd _eigenvalues;
  Eigen::MatrixXd _eigenvectors;

  Eigen::ArrayXi sort_index(Eigen::VectorXd &V) const;
  Eigen::MatrixXd get_initial_eigenvectors(Eigen::VectorXd &D, int size) const;

  Eigen::MatrixXd QR_ortho(Eigen::MatrixXd &A) const;
  Eigen::MatrixXd gramschmidt_ortho(Eigen::MatrixXd &A, int nstart) const;
  Eigen::VectorXd dpr_correction(Eigen::VectorXd &w, Eigen::VectorXd &A0,
                                 double lambda) const;
  Eigen::VectorXd olsen_correction(Eigen::VectorXd &r, Eigen::VectorXd &x,
                                   Eigen::VectorXd &D, double lambda) const;

  template <class MatrixReplacement>
  void update_projected_matrix(Eigen::MatrixXd &T, MatrixReplacement &A,
                               Eigen::MatrixXd &V) const {

    int size = V.rows();

    int old_dim = T.cols();
    int new_dim = V.cols();
    int nvec = new_dim - old_dim;

    T.conservativeResize(new_dim, new_dim);

    Eigen::MatrixXd _tmp = A * V.block(0, old_dim, size, nvec);

    T.block(0, old_dim, new_dim, nvec) = V.transpose() * _tmp;

    T.block(old_dim, 0, nvec, old_dim) =
        T.block(0, old_dim, old_dim, nvec).transpose();

    return;
  }
};
}  // namespace xtp
}  // namespace votca

#endif  // __VOTCA_TOOLS_DAVIDSON_SOLVER_H
