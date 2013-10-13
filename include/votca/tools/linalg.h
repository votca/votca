/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include "votca_config.h"

namespace votca { namespace tools {
    namespace ub = boost::numeric::ublas;

    /**
     * \brief solves A*x=b
     * @param x storage for x
     * @param A symmetric positive definite matrix for linear system
     * @param b inhomogenity
     * @param if A is not sysmetrix positive definite throws error code GSL_EDOM
     *
     * This function wrapps the cholesky linear system solver
     */
    void linalg_cholesky_solve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b);
    /**
     * \brief solves A*x=b
     * @param x storage for x
     * @param A matrix for linear equation system
     * @param b inhomogenity
     * @param residual if non-zero, residual will be stored here
     *
     * This function wrapps the qrsolver
     */
    void linalg_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::vector<double> *residual=NULL);

    /**
     * \brief solves A*x=b under the constraint B*x = 0
     * @param x storage for x
     * @param A matrix for linear equation system
     * @param b inhomogenity
     * @param constr constrained condition B (or is it the transposed one? check that)
     *
     * This function wrapps the qrsolver under constraints
     */
    void linalg_constrained_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::matrix<double> &constr);

    /**
     * \brief eigenvalues of a symmetric matrix A*x=E*x
     * @param A symmetric matrix 
     * @param E vector of eiganvalues
     * @param V matrix of eigenvalues
     * 
     * This function wrapps gsl_eigen_symmv
     * note that the eigenvalues/eigenvectors are UNSORTED 
     * 
     */
    bool linalg_eigenvalues_symmetric( ub::symmetric_matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V );
    
}}


#endif	/* __VOTCA_TOOLS_LINALG_H */

