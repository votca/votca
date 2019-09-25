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
#ifndef __VOTCA_TOOLS_LANCZOS_SOLVER_H
#define __VOTCA_TOOLS_LANCZOS_SOLVER_H

#include <chrono>
#include <iostream>
#include <stdexcept>

#include <boost/format.hpp>
#include <votca/xtp/eigen.h>
#include <votca/xtp/shiftinvert_operator.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/matrixfreeoperator.h>

#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/GenEigsSolver.h>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

/**
* \brief Use Lanczos algorithm to solve A*V=E*V

**/

class LanczosSolver {

 public:
  LanczosSolver(Logger &log);

  Eigen::ComputationInfo info() const { return _info; }

  Eigen::VectorXd eigenvalues() const { return this->_eigenvalues; }
  Eigen::MatrixXd eigenvectors() const { return this->_eigenvectors; }


  Eigen::ArrayXi argsort(const Eigen::VectorXd &V) const {
    /* \brief return the index of the sorted vector */
    Eigen::ArrayXi idx = Eigen::ArrayXi::LinSpaced(V.rows(), 0, V.rows() - 1);
    std::sort(idx.data(), idx.data() + idx.size(),
              [&](int i1, int i2) { return V[i1] < V[i2]; });
    return idx;
  }

  Eigen::MatrixXd sort_eigenvectors(const Eigen::MatrixXd &evec, Eigen::ArrayXi &idx)
  {
    Eigen::MatrixXd outvec = Eigen::MatrixXd::Zero(evec.rows(),evec.cols());
    for (int i=0; i<evec.cols(); i++){
      outvec.col(i) = evec.col(idx(i));
    }
    return outvec;
  }

  void set_iter_max(int N) { this->_iter_max = N; }
    
  void set_tolerance(std::string tol) {
    if (tol == "loose")
      this->_tol = 1E-4;
    else if (tol == "normal")
      this->_tol = 1E-6;
    else if (tol == "strict")
      this->_tol = 1E-8;
    else
      throw std::runtime_error(tol + " is not a valid Lanczos tolerance");
  }

  template <typename MatrixReplacement>
  void solve(const MatrixReplacement &A, int neigen) {

    std::chrono::time_point<std::chrono::system_clock> start =
        std::chrono::system_clock::now();
    int op_size = A.rows();
    PrintOptions(op_size);

    // declare the shift invert op
    ShiftInvertOperator<MatrixReplacement> sinv_op(A);

    // convergence criteria
    Eigen::Index ncv = 20;
    Eigen::Index nev = neigen;

    // simga
    double sigma = 0.;

    // solver
    Spectra::GenEigsRealShiftSolver<double, Spectra::LARGEST_REAL, 
      ShiftInvertOperator<MatrixReplacement>> eigs(&sinv_op, nev, ncv, sigma);

    //solve
    eigs.init();
    int nconv = eigs.compute(_iter_max, _tol);

    if (eigs.info() == Spectra::SUCCESSFUL)
    {
      this->_eigenvalues = eigs.eigenvalues().real();

      if (nev > 1) {
        std::sort(_eigenvalues.data(),_eigenvalues.data()+_eigenvalues.size());
        Eigen::ArrayXi _tmp_idx = LanczosSolver::argsort(_eigenvalues);
        Eigen::MatrixXd _tmp_evec = eigs.eigenvectors().real();
        this->_eigenvectors = LanczosSolver::sort_eigenvectors(_tmp_evec,_tmp_idx);
      }
       XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() <<"\nLanczos diagonalization converged in " << eigs.num_iterations()
         << flush;
    } else {
      this->_eigenvalues = Eigen::VectorXd::Zero(neigen);
      this->_eigenvectors = Eigen::MatrixXd::Zero(op_size,neigen);
      XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() <<"\nLanczos diagonalization failed" << flush;
    }
    _num_iterations = eigs.num_iterations();
    PrintTiming(start);
  }

  int num_iterations() {return _num_iterations;}

 private:
  Logger &_log;
  
  Eigen::VectorXd _eigenvalues;
  Eigen::MatrixXd _eigenvectors;
  double _tol = 1E-6;
  int _iter_max = 1000;
  int _num_iterations=0;

  Eigen::ComputationInfo _info = Eigen::ComputationInfo::NoConvergence;
  
  void PrintOptions(int op_size) const;
  void PrintTiming(
      const std::chrono::time_point<std::chrono::system_clock> &start) const;
  
};
}  // namespace xtp
}  // namespace votca

#endif  // __VOTCA_TOOLS_LANCZOS_SOLVER_H
