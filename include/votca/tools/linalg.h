/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_TOOLS_LINALG_H
#define	__VOTCA_TOOLS_LINALG_H
#include <votca/tools/eigen.h>

namespace votca { namespace tools {
 
    /**
     * \brief solves A*x=b under the constraint B*x = 0
     * @param x storage for x
     * @param A matrix for linear equation system
     * @param b inhomogenity
     * @param constr constrained condition
     *
     * This function implements the qrsolver under constraints
     */
    void linalg_constrained_qrsolve(Eigen::VectorXd &x, Eigen::MatrixXd &A,
            const Eigen::VectorXd &b, const Eigen::MatrixXd &constr);
  
  /**
     * \brief solves A*V=E*V for the first n eigenvalues
     * @param A symmetric matrix to diagonalize, is destroyed during iteration
     * @param E, eigenvalues
     * @param V, eigenvectors, each column is one eigenvector
     * @param nmax number of eigenvalues to return
     *
     * This function is only useful if MKL is used, wraps LAPACKE_ssyevx/LAPACKE_dsyevx
     */
    bool linalg_eigenvalues(Eigen::MatrixXd&A, Eigen::VectorXd &E, Eigen::MatrixXd&V , int nmax );
    bool linalg_eigenvalues(Eigen::MatrixXf&A, Eigen::VectorXf &E, Eigen::MatrixXf&V , int nmax );
   
   
}}



#endif	// __VOTCA_TOOLS_LINALG_H 

