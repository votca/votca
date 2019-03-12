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

#include <iostream>
#include <stdexcept>
#include <chrono>

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
  void set_tolerance(double eps) { this->tol = eps; }
  void set_max_search_space(int N) { this->max_search_space = N; }

  void set_correction(std::string method);

  Eigen::VectorXd eigenvalues() const { return this->_eigenvalues; }
  Eigen::MatrixXd eigenvectors() const { return this->_eigenvectors; }

  template <typename MatrixReplacement>
  void solve(MatrixReplacement &A, int neigen, int size_initial_guess = 0) {

    CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Davidson Solver"
                                 << flush;

    switch (this->davidson_correction) {

          case CORR::DPR:
            CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() 
            << " DPR Correction" << flush;
            break;

          case CORR::OLSEN:
            CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() 
            << " Olsen Correction" << flush;
            break;
        }
    CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() 
    << " Tolerance : " << this->tol << flush;

    double res_norm;
    double conv;
    int size = A.rows();
    bool has_converged = false;

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

    // initialize the guess eigenvector
    Eigen::VectorXd Adiag = A.diagonal();

    // target the lowest diagonal element
    Eigen::MatrixXd V =
    DavidsonSolver::_get_initial_eigenvectors(Adiag, size_initial_guess);

    // use a simple identity matrix
    //Eigen::MatrixXd V = Eigen::MatrixXd::Identity(size,size_initial_guess);

    // order the diagonal element
    // it seems to be needed or sometimes we
    // don't get the lowest eigenvalues .... t
    // std::sort(Adiag.data(),Adiag.data()+Adiag.size());

    Eigen::VectorXd lambda;  // eigenvalues hodlers
    Eigen::VectorXd old_val = Eigen::VectorXd::Zero(neigen);

    // temp varialbes
    Eigen::MatrixXd T, U, q;
    Eigen::VectorXd w, tmp;

    // project the matrix on the trial subspace
    T = A * V;
    T = V.transpose() * T;
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " iter\tSearch Space\tNorm" << flush;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_time;
    

    for (int iiter = 0; iiter < iter_max; iiter++) {

      start = std::chrono::system_clock::now();

      // diagonalize the small subspace
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
      lambda = es.eigenvalues();
      U = es.eigenvectors();

      // Ritz eigenvectors
      q = V.block(0, 0, V.rows(), search_space) * U;
      res_norm = 0.0;
      // correction vectors
      for (int j = 0; j < neigen; j++) {

        // residue vector
        w = A * q.col(j) - lambda(j) * q.col(j);
        res_norm += w.norm() / neigen;
        switch (this->davidson_correction) {

          case CORR::DPR:
            w = DavidsonSolver::_dpr_correction(w, Adiag, lambda(j));
            break;

          case CORR::OLSEN:
            tmp = q.col(j);
            w = DavidsonSolver::_olsen_correction(w, tmp, Adiag, lambda(j));
            break;
        }

        // append the correction vector to the search space
        V.conservativeResize(Eigen::NoChange, V.cols() + 1);
        V.col(V.cols() - 1) = w.normalized();
      }

      // Get the convergence criteria on the eigenvalues
      //conv = (lambda.head(neigen) - old_val).norm();

      // Print iteration data
      CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp()
                                   << format(" %1$4d \t %2$12d \t %3$4.2e") %
                                          iiter % search_space % res_norm
                                   << flush;
      // update
      search_space = V.cols();
      old_val = lambda.head(neigen);
      
      // break if converged
      if (res_norm < tol) {
        has_converged = true;
        break;
      }

      // check if we need to restart
      if (search_space > max_search_space or search_space > size) {
        V = q.block(0, 0, V.rows(), neigen);
        for (int j = 0; j < neigen; j++) {
          V.col(j).normalize();
        }
        search_space = neigen;

        // recompute the projected matrix
        T = A * V;
        T = V.transpose() * T;
      }

      // continue otherwise
      else {
        
        // orthogonalize the V vectors
        V = DavidsonSolver::_QR(V);

        // update the T matrix : avoid recomputing V.T A V
        // just recompute the element relative to the new eigenvectors
        DavidsonSolver::_update_projected_matrix<MatrixReplacement>(T, A, V);
      }
    }

    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << "-----------------------------------" << flush;
    if (!has_converged) {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << "- Warning : Davidson didn't converge ! "
          << flush;
      this->_eigenvalues = Eigen::VectorXd::Zero(neigen);
      this->_eigenvectors = Eigen::MatrixXd::Zero(size, neigen);
    } else {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << "- Davidson converged " << flush;

      // store the eigenvalues/eigenvectors
      this->_eigenvalues = lambda.head(neigen);
      this->_eigenvectors = q.block(0, 0, q.rows(), neigen);

      // normalize the eigenvectors
      for (int i = 0; i < neigen; i++) {
        this->_eigenvectors.col(i).normalize();
      }
    }
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << "-----------------------------------" << flush;
  }

 private:
  ctp::Logger &_log;
  int iter_max = 1000;
  double tol = 1E-3;
  int max_search_space = 100;

  enum CORR { DPR, OLSEN };
  CORR davidson_correction = CORR::DPR;

  Eigen::VectorXd _eigenvalues;
  Eigen::MatrixXd _eigenvectors;

  Eigen::ArrayXi _sort_index(Eigen::VectorXd &V) const;
  Eigen::MatrixXd _get_initial_eigenvectors(Eigen::VectorXd &D, int size) const;

  Eigen::MatrixXd _QR(Eigen::MatrixXd &A) const;

  Eigen::VectorXd _dpr_correction(Eigen::VectorXd &w, Eigen::VectorXd &A0,
                                  double lambda) const;
  Eigen::VectorXd _olsen_correction(Eigen::VectorXd &r, Eigen::VectorXd &x,
                                    Eigen::VectorXd &D, double lambda) const;

  template <class MatrixReplacement>
  void _update_projected_matrix(Eigen::MatrixXd &T, MatrixReplacement &A,
                                Eigen::MatrixXd &V) const {
    int nvec_old = T.cols();
    int nvec = V.cols();
    int nnew_vec = nvec - nvec_old;

    Eigen::MatrixXd _tmp = A * V.block(0, nvec_old, V.rows(), nnew_vec);
    T.conservativeResize(nvec, nvec);
    T.block(0, nvec_old, nvec, nnew_vec) = V.transpose() * _tmp;
    T.block(nvec_old, 0, nnew_vec, nvec_old) =
        T.block(0, nvec_old, nvec_old, nnew_vec).transpose();
    return;
  }
};
}  // namespace xtp
}  // namespace votca

#endif  // __VOTCA_TOOLS_DAVIDSON_SOLVER_H
