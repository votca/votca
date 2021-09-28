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

#pragma once
#ifndef VOTCA_XTP_DAVIDSONSOLVER_H
#define VOTCA_XTP_DAVIDSONSOLVER_H

// Standard includes
#include <chrono>
#include <iostream>
#include <stdexcept>

// Third party includes
#include <boost/format.hpp>

// Local VOTCA includes
#include "eigen.h"
#include "logger.h"

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

/**
 * \brief Use Davidson algorithm to solve A*V=E*V
 *
 * For a rough overview https://joshuagoings.com/2013/08/23/davidsons-method/
 * For the symmetric case we simply extract the smallest eigenvalues
 * For the non-symmetric case we need the smallest non-negative eigenvalues
 * These are harder to extract, because iterative methods tend to converge
 *towards extreme eigenvalues, thus we use harmonic ritz values.
 **/

class DavidsonSolver {

 public:
  DavidsonSolver(Logger &log);

  void set_iter_max(Index N) { this->_iter_max = N; }
  void set_max_search_space(Index N) { this->_max_search_space = N; }
  void set_tolerance(std::string tol);
  void set_correction(std::string method);
  void set_size_update(std::string update_size);
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
    _restart_size = size_initial_guess;

    // get the diagonal of the operator
    this->_Adiag = A.diagonal();

    // target the lowest diagonal element
    ProjectedSpace proj = initProjectedSpace(neigen, size_initial_guess);
    RitzEigenPair rep;
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " iter\tSearch Space\tNorm" << flush;

    for (_i_iter = 0; _i_iter < _iter_max; _i_iter++) {

      updateProjection(A, proj);

      rep = getRitzEigenPairs(proj);

      bool converged = checkConvergence(rep, proj, neigen);

      printIterationData(rep, proj, neigen);

      bool last_iter = _i_iter == (_iter_max - 1);

      if (converged) {
        storeConvergedData(rep, neigen);
        break;
      } else if (last_iter) {
        storeNotConvergedData(rep, proj.root_converged, neigen);
        break;
      }
      Index extension_size = extendProjection(rep, proj);
      bool do_restart = (proj.search_space() > _max_search_space);

      if (do_restart) {
        restart(rep, proj, extension_size);
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
  Index _restart_size = 0;
  enum CORR { DPR, OLSEN };
  CORR _davidson_correction = CORR::DPR;

  enum UPDATE { MIN, SAFE, MAX };
  UPDATE _davidson_update = UPDATE::SAFE;

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

    // These are only used for harmonic ritz in the non-hermitian case
    Eigen::MatrixXd AAV;  // A*A*V
    Eigen::MatrixXd B;    // V.T *A*A*V
  };

  template <typename MatrixReplacement>
  void updateProjection(const MatrixReplacement &A,
                        ProjectedSpace &proj) const {

    if (_i_iter == 0) {
      proj.AV = A * proj.V;
      proj.T = proj.V.transpose() * proj.AV;
      if (_matrix_type == MATRIX_TYPE::HAM) {
        proj.AAV = A * proj.AV;
        proj.B = proj.V.transpose() * proj.AAV;
      }

    } else {
      /* if we use a Gram Schmid(GS) orthogonalisation we do not have to
      recompute the entire projection as GS doesn't change the original
      subspace*/
      Index old_dim = proj.AV.cols();
      Index new_dim = proj.V.cols();
      Index nvec = new_dim - old_dim;
      proj.AV.conservativeResize(Eigen::NoChange, new_dim);
      proj.AV.rightCols(nvec) = A * proj.V.rightCols(nvec);

      proj.T.conservativeResize(new_dim, new_dim);
      proj.T.rightCols(nvec) = proj.V.transpose() * proj.AV.rightCols(nvec);

      if (_matrix_type == MATRIX_TYPE::SYMM) {
        proj.T.bottomLeftCorner(nvec, old_dim) =
            proj.T.topRightCorner(old_dim, nvec).transpose();

      } else {
        proj.T.bottomLeftCorner(nvec, old_dim) =
            proj.V.rightCols(nvec).transpose() * proj.AV.leftCols(old_dim);

        proj.AAV.conservativeResize(Eigen::NoChange, new_dim);
        proj.AAV.rightCols(nvec) = A * proj.AV.rightCols(nvec);
        proj.B.conservativeResize(new_dim, new_dim);
        proj.B.rightCols(nvec) = proj.V.transpose() * proj.AAV.rightCols(nvec);
        proj.B.bottomLeftCorner(nvec, old_dim) =
            proj.V.rightCols(nvec).transpose() * proj.AAV.leftCols(old_dim);
      }
    }
  }
  RitzEigenPair getRitzEigenPairs(const ProjectedSpace &proj) const;

  Eigen::MatrixXd qr(const Eigen::MatrixXd &A) const;
  RitzEigenPair getHarmonicRitz(const ProjectedSpace &proj) const;

  RitzEigenPair getRitz(const ProjectedSpace &proj) const;

  Index getSizeUpdate(Index neigen) const;

  void checkOptions(Index operator_size);

  void printOptions(Index operator_size) const;

  void printTiming(
      const std::chrono::time_point<std::chrono::system_clock> &start) const;

  void printIterationData(const RitzEigenPair &rep, const ProjectedSpace &proj,
                          Index neigen) const;

  ArrayXl argsort(const Eigen::VectorXd &V) const;

  Eigen::MatrixXd setupInitialEigenvectors(Index size_initial_guess) const;

  Eigen::MatrixXd extract_vectors(const Eigen::MatrixXd &V,
                                  const ArrayXl &idx) const;

  void orthogonalize(Eigen::MatrixXd &V, Index nupdate) const;
  void gramschmidt(Eigen::MatrixXd &A, Index nstart) const;

  Eigen::VectorXd computeCorrectionVector(const Eigen::VectorXd &qj,
                                          double lambdaj,
                                          const Eigen::VectorXd &Aqj) const;
  Eigen::VectorXd dpr(const Eigen::VectorXd &r, double lambda) const;
  Eigen::VectorXd olsen(const Eigen::VectorXd &r, const Eigen::VectorXd &x,
                        double lambda) const;

  ProjectedSpace initProjectedSpace(Index neigen,
                                    Index size_initial_guess) const;

  Index extendProjection(const RitzEigenPair &rep, ProjectedSpace &proj);

  bool checkConvergence(const RitzEigenPair &rep, ProjectedSpace &proj,
                        Index neigen);

  void restart(const RitzEigenPair &rep, ProjectedSpace &proj,
               Index newtestvectors) const;

  void storeConvergedData(const RitzEigenPair &rep, Index neigen);

  void storeNotConvergedData(const RitzEigenPair &rep,
                             std::vector<bool> &root_converged, Index neigen);

  void storeEigenPairs(const RitzEigenPair &rep, Index neigen);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DAVIDSONSOLVER_H
