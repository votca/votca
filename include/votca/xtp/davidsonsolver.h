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

#pragma once
#ifndef __VOTCA_TOOLS_DAVIDSON_SOLVER_H
#define __VOTCA_TOOLS_DAVIDSON_SOLVER_H

#include <chrono>
#include <iostream>
#include <stdexcept>

#include <boost/format.hpp>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

/**
* \brief Use Davidson algorithm to solve A*V=E*V

**/

class DavidsonSolver {

 public:
  DavidsonSolver(Logger &log);

  void set_iter_max(int N) { this->_iter_max = N; }

  void set_max_search_space(int N) { this->_max_search_space = N; }

  void set_tolerance(std::string tol);
  void set_correction(std::string method);
  void set_ortho(std::string method);
  void set_size_update(std::string method);
  void set_matrix_type(std::string mt);

  Eigen::ComputationInfo info() const { return _info; }

  Eigen::VectorXd eigenvalues() const { return this->_eigenvalues; }
  Eigen::MatrixXd eigenvectors() const { return this->_eigenvectors; }
  Eigen::MatrixXd residues() const { return this->_res; }

  template <typename MatrixReplacement>
  void solve(const MatrixReplacement &A, int neigen,
             int size_initial_guess = 0) {

    std::chrono::time_point<std::chrono::system_clock> start =
        std::chrono::system_clock::now();
    int op_size = A.rows();
    PrintOptions(op_size);
    //. search space exceeding the system size
    if (_max_search_space > op_size) {
      XTP_LOG_SAVE(logDEBUG, _log)
          << TimeStamp() << " Warning Max search space (" << _max_search_space
          << ") larger than system size (" << op_size << ")" << flush;
      _max_search_space = op_size;
      XTP_LOG_SAVE(logDEBUG, _log)
          << TimeStamp() << " Max search space set to " << op_size << flush;
    }

    // initial guess size
    if (size_initial_guess == 0) {
      size_initial_guess = 2 * neigen;
    }

    int search_space = size_initial_guess;
    int size_restart = size_initial_guess;
    int size_update = get_size_update(neigen);

    std::vector<bool> root_converged = std::vector<bool>(size_update, false);

    // initialize the guess eigenvector
    Eigen::VectorXd Adiag = A.diagonal();

    // target the lowest diagonal element
    Eigen::MatrixXd V = SetupInitialEigenvectors(Adiag, size_initial_guess);
    Eigen::MatrixXd q;
    Eigen::MatrixXd U;
    Eigen::VectorXd lambda;
    Eigen::MatrixXd AV;
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " iter\tSearch Space\tNorm" << flush;

    // Start of the main iteration loop
    for (int iiter = 0; iiter < _iter_max; iiter++) {

      Eigen::MatrixXd T;

      // check if we need to restart
      bool restart_required =
          search_space > _max_search_space || search_space > op_size;
      if (iiter == 0 || _davidson_ortho == ORTHO::QR) {
        AV = A * V;
        T = V.transpose() * AV;
      } else if (restart_required) {
        V = q.leftCols(size_restart);
        V.colwise().normalize();
        AV = AV * U.leftCols(size_restart);  // corresponds to replacing V with
                                             // q.leftCols
        search_space = size_restart;
        T = V.transpose() * AV;
      } else if (_davidson_ortho == ORTHO::GS) {  // GS orthogonalisation leaves
                                                  // old eigenvectors untouched,
                                                  // so we only have to update
                                                  // the new ones
        int size = V.rows();
        int old_dim = T.cols();
        int new_dim = V.cols();
        int nvec = new_dim - old_dim;
        AV.conservativeResize(Eigen::NoChange, new_dim);
        AV.block(0, old_dim, size, nvec) = A * V.block(0, old_dim, size, nvec);
        Eigen::MatrixXd VAV = V.transpose() * AV.block(0, old_dim, size, nvec);
        T.conservativeResize(new_dim, new_dim);
        T.block(0, old_dim, new_dim, nvec) = VAV;
        T.block(old_dim, 0, nvec, old_dim) = VAV.topRows(old_dim).transpose();
      }

      // diagonalize the small subspace
      switch (this->_matrix_type) {
        case TYPE::SYMM: {
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
          lambda = es.eigenvalues();
          U = es.eigenvectors();
          break;
        }
        case TYPE::HAM: {

          Eigen::EigenSolver<Eigen::MatrixXd> es(T);
          lambda = es.eigenvalues().real();
          U = es.eigenvectors().real();

          Eigen::ArrayXi reorder_idx = 
            index_window(lambda,size_update,0.0,0.25);
          lambda = reorder_idx.unaryExpr(lambda);
          U = extract_eigenvectors(U, reorder_idx);  
          U.colwise().normalize();
          break;  
        }
      }
      
    
      q = V * U;  // Ritz vectors
      Eigen::MatrixXd r =
          AV * U - q * lambda.asDiagonal();  // compute residual=A*Q - lambda Q

      Eigen::ArrayXd res_norm = Eigen::ArrayXd::Zero(size_update);
      int nupdate = 0;

      for (int j = 0; j < size_update; j++) {
        // skip the root that have already converged
        if (root_converged[j]) {
          continue;
        }
        nupdate++;

        res_norm[j] = r.col(j).norm();
        
        // residue vector
        Eigen::VectorXd w =
            ComputeCorrectionVector(Adiag, q.col(j), lambda(j), r.col(j));

        // append the correction vector to the search space
        V.conservativeResize(Eigen::NoChange, V.cols() + 1);
        V.rightCols<1>() = w.normalized();

        // track converged root
        root_converged[j] = (res_norm[j] < _tol);
      }

      // Print iteration data
      int converged_roots = 0;
      for (int i = 0; i < neigen; i++) {
        converged_roots += root_converged[i];
      }
      double percent_converged = 100 * double(converged_roots) / double(neigen);
      XTP_LOG_SAVE(logDEBUG, _log)
          << TimeStamp()
          << format(" %1$4d %2$12d \t %3$4.2e \t %4$5.2f%% converged") % iiter %
                 search_space % res_norm.head(neigen).maxCoeff() %
                 percent_converged
          << flush;

      // update
      search_space = V.cols();

      // coverged
      bool converged;
      switch (this->_matrix_type){
        case TYPE::SYMM: {
           converged = (res_norm.head(neigen) < _tol).all();
        }
        case TYPE::HAM: {
          Eigen::ArrayXi idx = index_window(lambda,neigen,0,0);
          converged = (idx.unaryExpr(res_norm) < _tol).all();
        }
      }

      if (converged) {
        _info = Eigen::ComputationInfo::Success;
      }

      bool last_iter = iiter == (_iter_max - 1);

      // break if converged or last
      if (converged || last_iter) {

        // store the eigenvalues/eigenvectors
        Eigen::ArrayXi idx = index_window(lambda,neigen,0,0);
        this->_eigenvalues = idx.unaryExpr(lambda);
        this->_eigenvectors = extract_eigenvectors(q, idx);
        this->_res = idx.unaryExpr(res_norm);

        this->_eigenvectors.colwise().normalize();
        if (last_iter && !converged) {
          XTP_LOG_SAVE(logDEBUG, _log)
              << TimeStamp() << "- Warning : Davidson " << percent_converged
              << "% converged after " << _iter_max << " iterations." << flush;
          _info = Eigen::ComputationInfo::NoConvergence;
          for (int i = 0; i < neigen; i++) {
            if (!root_converged[i]) {
              _eigenvalues(i) = 0;
              _eigenvectors.col(i) = Eigen::VectorXd::Zero(op_size);
            }
          }
        }
        break;
      }

      switch (this->_davidson_ortho) {
        case ORTHO::GS:
          V = DavidsonSolver::gramschmidt_ortho(V, V.cols() - nupdate);
          break;
        case ORTHO::QR:
          V = DavidsonSolver::QR_ortho(V);
          break;
      }
    }

    PrintTiming(start);
  }

 private:
  Logger &_log;
  int _iter_max = 50;
  double _tol = 1E-4;
  int _max_search_space = 1000;

  enum CORR { DPR, OLSEN };
  CORR _davidson_correction = CORR::DPR;

  enum UPDATE { MIN, SAFE, MAX };
  UPDATE _davidson_update = UPDATE::SAFE;

  enum ORTHO { GS, QR };
  ORTHO _davidson_ortho = ORTHO::GS;

  enum TYPE { SYMM, HAM };  
  TYPE _matrix_type = TYPE::SYMM;

  Eigen::VectorXd _eigenvalues;
  Eigen::MatrixXd _eigenvectors;
  Eigen::VectorXd _res;

  Eigen::ComputationInfo _info = Eigen::ComputationInfo::NoConvergence;
  int get_size_update(int neigen) const;

  void PrintOptions(int op_size) const;
  void PrintTiming(
      const std::chrono::time_point<std::chrono::system_clock> &start) const;
  Eigen::VectorXd ComputeCorrectionVector(const Eigen::VectorXd &Adiag,
                                          const Eigen::VectorXd &qj,
                                          double lambdaj,
                                          const Eigen::VectorXd &Aqj) const;

  Eigen::ArrayXi argsort(const Eigen::VectorXd &V) const;
  Eigen::MatrixXd SetupInitialEigenvectors(Eigen::VectorXd &D, int size) const;

  Eigen::ArrayXi index_window(const Eigen::VectorXd &V, 
    int size_update, double target_min_val, double perc_below) const;

  Eigen::MatrixXd extract_eigenvectors(const Eigen::MatrixXd &V, 
    const Eigen::ArrayXi &idx) const;

  Eigen::MatrixXd QR_ortho(const Eigen::MatrixXd &A) const;
  Eigen::MatrixXd gramschmidt_ortho(const Eigen::MatrixXd &A, int nstart);
  Eigen::VectorXd dpr_correction(const Eigen::VectorXd &w,
                                 const Eigen::VectorXd &A0,
                                 double lambda) const;
  Eigen::VectorXd olsen_correction(const Eigen::VectorXd &r,
                                   const Eigen::VectorXd &x,
                                   const Eigen::VectorXd &D,
                                   double lambda) const;
};
}  // namespace xtp
}  // namespace votca

#endif  // __VOTCA_TOOLS_DAVIDSON_SOLVER_H
