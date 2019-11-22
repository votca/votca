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

  void set_iter_max(Index N) { this->_iter_max = N; }
  void set_max_search_space(Index N) { this->_max_search_space = N; }
  void set_tolerance(std::string tol);
  void set_correction(std::string method);
  void set_ortho(std::string method);
  void set_size_update(std::string method);
  void set_matrix_type(std::string mt);

  Eigen::ComputationInfo info() const { return _info; }
  Eigen::VectorXd eigenvalues() const { return this->_eigenvalues; }
  Eigen::MatrixXd eigenvectors() const { return this->_eigenvectors; }
  Eigen::MatrixXd residues() const { return this->_res; }
  Index num_iterations() const { return this->_i_iter; }

  template <typename MatrixReplacement>
  void solve(const MatrixReplacement &A, Index neigen,
             Index size_initial_guess = 0) {

    if (_max_search_space < neigen) {
      _max_search_space = neigen * 5;
    }
    std::chrono::time_point<std::chrono::system_clock> start =
        std::chrono::system_clock::now();
    Index op_size = A.rows();

    checkOptions(op_size);
    printOptions(op_size);

    // initial guess size
    if (size_initial_guess == 0) {
      size_initial_guess = 2 * neigen;
    }

    // get the diagonal of the operator
    this->_Adiag = A.diagonal();

    // target the lowest diagonal element
    ProjectedSpace proj = initProjectedSpace(neigen, size_initial_guess);
    RitzEigenPair rep;

    for (_i_iter = 0; _i_iter < _iter_max; _i_iter++) {

      bool do_restart = (proj.search_space() > _max_search_space);

      if (do_restart) {
        restart(rep, proj, size_initial_guess);

      } else {
        updateProjection(A, proj);
        rep = getGuessVectors(A, proj);
        extendProjection(rep, proj);
      }

      printIterationData(rep, proj, neigen);

      bool converged = (rep.res_norm().head(neigen) < _tol).all();
      bool last_iter = _i_iter == (_iter_max - 1);

      if (converged) {
        storeConvergedData(rep, neigen);
        break;
      } else if (last_iter) {
        storeNotConvergedData(rep, proj.root_converged, neigen);
        break;
      }
    }

    printTiming(start);
  }

 private:
  Logger &_log;
  Index _iter_max = 50;
  Index _i_iter = 0;
  double _tol = 1E-4;
  Index _max_search_space = 0;
  Eigen::VectorXd _Adiag;

  enum CORR { DPR, OLSEN };
  CORR _davidson_correction = CORR::DPR;

  enum UPDATE { MIN, SAFE, MAX };
  UPDATE _davidson_update = UPDATE::SAFE;

  enum ORTHO { GS, QR };
  ORTHO _davidson_ortho = ORTHO::GS;

  enum MATRIX_TYPE { SYMM, HAM };
  MATRIX_TYPE _matrix_type = MATRIX_TYPE::SYMM;

  Eigen::VectorXd _eigenvalues;
  Eigen::MatrixXd _eigenvectors;
  Eigen::VectorXd _res;
  Eigen::ComputationInfo _info = Eigen::ComputationInfo::NoConvergence;

  struct RitzEigenPair {
    Eigen::VectorXd lambda;  // eigenvalues
    Eigen::MatrixXd q;       // Ritz (or harmonic Ritz) eigenvectors
    Eigen::MatrixXd U;       // eigenvectors of the small subspace
    Eigen::MatrixXd res;     // residues of the pairs
    Eigen::ArrayXd res_norm() const {
      return res.colwise().norm();
    }  // norm of the residues
  };

  struct ProjectedSpace {
    Eigen::MatrixXd V;   // basis of vectors
    Eigen::MatrixXd AV;  // A * V
    Eigen::MatrixXd T;   // V.T * A * V
    Index search_space() const {
      return V.cols();
    };                  // size of the projection i.e. number of cols in V
    Index size_update;  // size update ...
    std::vector<bool> root_converged;  // keep track of which root have onverged
  };

  template <typename MatrixReplacement>
  void updateProjection(const MatrixReplacement &A,
                        ProjectedSpace &proj) const {

    if (_i_iter == 0 || _davidson_ortho == ORTHO::QR) {
      /* if we use QR we need to recompute the entire projection
      since QR will modify original subspace*/
      proj.AV = A * proj.V;
      proj.T = proj.V.transpose() * proj.AV;

    } else if (_davidson_ortho == ORTHO::GS) {
      /* if we use a GS ortho we do not have to recompute
      the entire projection as GS doesn't change the original subspace*/
      Index old_dim = proj.T.cols();
      Index new_dim = proj.V.cols();
      Index nvec = new_dim - old_dim;
      proj.AV.conservativeResize(Eigen::NoChange, new_dim);
      proj.AV.rightCols(nvec) = A * proj.V.rightCols(nvec);
      Eigen::MatrixXd VAV = proj.V.transpose() * proj.AV.rightCols(nvec);
      proj.T.conservativeResize(new_dim, new_dim);
      proj.T.rightCols(nvec) = VAV;
      proj.T.bottomLeftCorner(nvec, old_dim) =
          proj.T.topRightCorner(old_dim, nvec).transpose();
    }
  }

  template <typename MatrixReplacement>
  RitzEigenPair getGuessVectors(const MatrixReplacement &A,
                                const ProjectedSpace &proj) const {
    // get the ritz vectors
    switch (this->_matrix_type) {
      case MATRIX_TYPE::SYMM: {
        return getRitz(proj);
      }
      case MATRIX_TYPE::HAM: {
        return getHarmonicRitz(A, proj);
      }
    }
    return RitzEigenPair();
  }

  template <typename MatrixReplacement>
  RitzEigenPair getHarmonicRitz(const MatrixReplacement &A,
                                const ProjectedSpace &proj) const {

    /* Compute the Harmonic Ritz vector following
     * Computing Interior Eigenvalues of Large Matrices
     * Ronald B Morgan
     * LINEAR ALGEBRA AND ITS APPLICATIONS 154-156:289-309 (1991)
     * https://cpb-us-w2.wpmucdn.com/sites.baylor.edu/dist/e/71/files/2015/05/InterEvals-1vgdz91.pdf
     */

    RitzEigenPair rep;
    Eigen::MatrixXd B = proj.V.transpose() * (A * proj.AV);

    bool return_eigenvectors = true;
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges(proj.T, B,
                                                       return_eigenvectors);

    ArrayXl idx = DavidsonSolver::argsort(ges.eigenvalues().real());
    // smallest to largest
    idx = idx.reverse();
    // we need the largest values, because this is the inverse value, so reverse
    // list

    rep.U = DavidsonSolver::extract_vectors(ges.eigenvectors().real(), idx);
    rep.U.colwise().normalize();
    rep.lambda = (rep.U.transpose() * proj.T * rep.U).diagonal();

    rep.q = proj.V * rep.U;  // Ritz vectors
    rep.res = proj.AV * rep.U - rep.q * rep.lambda.asDiagonal();  // residues
    return rep;
  }

  RitzEigenPair getRitz(const ProjectedSpace &proj) const;

  Index getSizeUpdate(Index neigen) const;

  void checkOptions(Index operator_size);

  void printOptions(Index operator_size) const;

  void printTiming(
      const std::chrono::time_point<std::chrono::system_clock> &start) const;

  void printIterationData(const RitzEigenPair &rep, const ProjectedSpace &proj,
                          Index neigen) const;

  ArrayXl argsort(const Eigen::VectorXd &V) const;

  Eigen::MatrixXd setupInitialEigenvectors(Index size) const;

  Eigen::MatrixXd extract_vectors(const Eigen::MatrixXd &V,
                                  const ArrayXl &idx) const;

  Eigen::MatrixXd orthogonalize(const Eigen::MatrixXd &V, Index nupdate);
  Eigen::MatrixXd qr(const Eigen::MatrixXd &A) const;
  Eigen::MatrixXd gramschmidt(const Eigen::MatrixXd &A, Index nstart);

  Eigen::VectorXd computeCorrectionVector(const Eigen::VectorXd &qj,
                                          double lambdaj,
                                          const Eigen::VectorXd &Aqj) const;
  Eigen::VectorXd dpr(const Eigen::VectorXd &w, double lambda) const;
  Eigen::VectorXd olsen(const Eigen::VectorXd &r, const Eigen::VectorXd &x,
                        double lambda) const;

  ProjectedSpace initProjectedSpace(Index neigen,
                                    Index size_initial_guess) const;

  void extendProjection(const RitzEigenPair &rep, ProjectedSpace &proj);

  void restart(const RitzEigenPair &rep, ProjectedSpace &proj,
               Index size_restart) const;

  void storeConvergedData(const RitzEigenPair &rep, Index neigen);

  void storeNotConvergedData(const RitzEigenPair &rep,
                             std::vector<bool> &root_converged, Index neigen);

  void storeEigenPairs(const RitzEigenPair &rep, Index neigen);
};

}  // namespace xtp
}  // namespace votca

#endif  // __VOTCA_TOOLS_DAVIDSON_SOLVER_H
