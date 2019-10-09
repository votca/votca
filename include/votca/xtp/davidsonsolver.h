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
  int num_iterations() const { return this->_num_iter; }
  
  template <typename MatrixReplacement>
  void solve(const MatrixReplacement &A, int neigen,
             int size_initial_guess = 0) {

    std::chrono::time_point<std::chrono::system_clock> start =
        std::chrono::system_clock::now();
    int op_size = A.rows();

    checkOptions(op_size);
    printOptions(op_size);

    // initial guess size
    if (size_initial_guess == 0) {
      size_initial_guess = 2 * neigen;
    }

    int size_update = getSizeUpdate(neigen);
    std::vector<bool> root_converged = std::vector<bool>(size_update, false);

    // initialize the guess eigenvector
    Eigen::VectorXd Adiag = A.diagonal();

    // target the lowest diagonal element
    ProjectedSpace proj;
    proj.V = setupInitialEigenvectors(Adiag, size_initial_guess);
    proj.search_space = proj.V.cols();

    // ritz eigenpair
    RitzEigenPair rep;

    // Start of the main iteration loop
    int nupdate;
    for (int iiter = 0; iiter < _iter_max; iiter++) {
      
      // check if we need to restart
      bool restart_required = proj.search_space > _max_search_space;

      if (iiter == 0 || _davidson_ortho == ORTHO::QR) {
        proj.AV = A * proj.V;
        proj.T = proj.V.transpose() * proj.AV;
      } 
      else if (restart_required) {
        restart(rep,proj,size_initial_guess);

      } else if (_davidson_ortho == ORTHO::GS) {  
        updateProjection(A, proj);
      }

      // get the ritz vectors
      rep = computeRitzEigenPairs(proj,size_update);
      
      // etend the subspace
      nupdate = extendProjection(Adiag,rep,proj,root_converged,size_update);

      // Print iteration data
      printIterationData(root_converged, rep.res_norm, neigen, proj.search_space, iiter);

      // converged
      bool converged = checkConvergence(rep,neigen);
      bool last_iter = iiter == (_iter_max - 1);

      // break if converged or last
      if (converged) {
        storeConvergedData(rep, neigen, iiter);
        break;
      } else if (last_iter) {
        storeNotConvergedData(rep,root_converged, neigen);
        break;
      }

      proj.V = orthogonalize(proj.V,nupdate);
    }

    printTiming(start);
  }

 private:

  Logger &_log;
  int _iter_max = 50;
  int _num_iter = 0;
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

  struct RitzEigenPair {
    Eigen::VectorXd lambda; //eigenvalues
    Eigen::MatrixXd q; // Ritz (or harmonic Ritz) eigenvectors
    Eigen::MatrixXd U; //eigenvectors of the small subspace
    Eigen::MatrixXd res; // residues of the pairs
    Eigen::ArrayXd res_norm; // norm of the reisidues
  };

  struct ProjectedSpace {
    Eigen::MatrixXd V; // basis of vectors
    Eigen::MatrixXd AV; // A * V
    Eigen::MatrixXd T; // V.T * A * V
    int search_space; //size of the projection i.e. number of cols in V
  };

  template <typename MatrixReplacement>
  void updateProjection(const MatrixReplacement &A, ProjectedSpace &proj) const {
      /* if we use a GS ortho we do not have to recompute 
      the entire projection as GS doesn't change the original subspace*/
      int size = proj.V.rows();
      int old_dim = proj.T.cols();
      int new_dim = proj.V.cols();
      int nvec = new_dim - old_dim;
      proj.AV.conservativeResize(Eigen::NoChange, new_dim);
      proj.AV.block(0, old_dim, size, nvec) = A * proj.V.block(0, old_dim, size, nvec);
      Eigen::MatrixXd VAV = proj.V.transpose() * proj.AV.block(0, old_dim, size, nvec);
      proj.T.conservativeResize(new_dim, new_dim);
      proj.T.block(0, old_dim, new_dim, nvec) = VAV;
      proj.T.block(old_dim, 0, nvec, old_dim) = VAV.topRows(old_dim).transpose();
  }

  int getSizeUpdate(int neigen) const;

  void checkOptions(int op_size);

  void printOptions(int op_size) const;

  void printTiming(
      const std::chrono::time_point<std::chrono::system_clock> &start) const;

  void printIterationData(std::vector<bool> const &root_converged,
    Eigen::ArrayXd const &res, int neigen, int search_space, int iiter) const;


  Eigen::ArrayXi argsort(const Eigen::VectorXd &V) const;
  Eigen::MatrixXd setupInitialEigenvectors(Eigen::VectorXd &D, int size) const;
  RitzEigenPair computeRitzEigenPairs ( const ProjectedSpace &proj, int size_update);

  Eigen::ArrayXi index_window(const Eigen::VectorXd &V, 
    int size_update, double target_min_val, double perc_below) const;

  Eigen::MatrixXd extract_eigenvectors(const Eigen::MatrixXd &V, 
    const Eigen::ArrayXi &idx) const;

  Eigen::MatrixXd orthogonalize(const Eigen::MatrixXd &V, int nupdate);
  Eigen::MatrixXd qr(const Eigen::MatrixXd &A) const;
  Eigen::MatrixXd gs(const Eigen::MatrixXd &A, int nstart);

  Eigen::VectorXd computeCorrectionVector(const Eigen::VectorXd &Adiag,
                                          const Eigen::VectorXd &qj,
                                          double lambdaj,
                                          const Eigen::VectorXd &Aqj) const;
  Eigen::VectorXd dpr(const Eigen::VectorXd &w,
                                 const Eigen::VectorXd &A0,
                                 double lambda) const;
  Eigen::VectorXd olsen(const Eigen::VectorXd &r,
                                   const Eigen::VectorXd &x,
                                   const Eigen::VectorXd &D,
                                   double lambda) const;

  int extendProjection(Eigen::VectorXd &Adiag, RitzEigenPair &rep, ProjectedSpace &proj, 
    std::vector<bool> &root_converged, int size_update);

  void restart (const RitzEigenPair &rep, ProjectedSpace &proj, int size_restart) const;

  bool checkConvergence(const RitzEigenPair &rep, int neigen);

  void storeConvergedData(const RitzEigenPair &rep, int neigen, int iiter);

  void storeNotConvergedData(const RitzEigenPair &rep, std::vector<bool> &root_converged, 
                             int neigen);

  void storeEigenPairs(const RitzEigenPair &rep, int neigen);

};

}  // namespace xtp
}  // namespace votca

#endif  // __VOTCA_TOOLS_DAVIDSON_SOLVER_H
