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

  template <typename MatrixReplacement>
  void solve(const MatrixReplacement &A, int neigen) {

    std::chrono::time_point<std::chrono::system_clock> start =
        std::chrono::system_clock::now();

    // declare the shift invert op
    SINV_OPERATOR<MatrixReplacement> sinv_op(A);

    // convergence criteria
    Eigen::Index ncv = 2*neigen+1;
    Eigen::Index nev = neigen;

    // simga
    double sigma = 0.;

    // solver
    // Spectra::GenEigsRealShiftSolver<double, Spectra::SMALLEST_REAL, 
    //   SINV_OPERATOR<MatrixReplacement> > eigs(&sinv_op, nev, ncv, sigma);


    Spectra::GenEigsSolver<double, Spectra::SMALLEST_REAL, 
      SINV_OPERATOR<MatrixReplacement> > eigs(&sinv_op, nev, ncv);

    eigs.init();
    int nconv = eigs.compute();

    if (eigs.info() == Spectra::SUCCESSFUL)
    {
      this->_eigenvalues = eigs.eigenvalues().real();
      this->_eigenvectors = eigs.eigenvectors().real();
    }

    PrintTiming(start);
  }

 private:
  Logger &_log;
  
  Eigen::VectorXd _eigenvalues;
  Eigen::MatrixXd _eigenvectors;
  Eigen::ComputationInfo _info = Eigen::ComputationInfo::NoConvergence;
  
  void PrintOptions(int op_size) const;
  void PrintTiming(
      const std::chrono::time_point<std::chrono::system_clock> &start) const;
  
};
}  // namespace xtp
}  // namespace votca

#endif  // __VOTCA_TOOLS_LANCZOS_SOLVER_H
