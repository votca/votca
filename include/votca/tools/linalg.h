/* 
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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
#include <votca/tools/votca_config.h>
#if defined(GSL)
    #include "votca_gsl_boost_ublas_matrix_prod.h"
#elif defined(MKL)
    #include "mkl_boost_ublas_matrix_prod.hpp"
#endif

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

namespace votca { namespace tools {
    namespace ub = boost::numeric::ublas;

    /**
     * \brief inverts A
     * @param A symmetric positive definite matrix
     * @param V inverse matrix
     *
     * This function wraps the inversion of a matrix
     */
    void linalg_invert( ub::matrix<double> &A, ub::matrix<double> &V );
    void linalg_invert( ub::matrix<float> &A, ub::matrix<float> &V );
 
    /**
     * \brief determines Cholesky decomposition of matrix A
     * @param A symmetric positive definite matrix
     *
     * This function wraps the Cholesky decomposition
     */
    void linalg_cholesky_decompose( ub::matrix<double> &A );
    void linalg_cholesky_decompose( ub::matrix<float> &A );
    /**
     * \brief solves A*x=b
     * @param x storage for x
     * @param A symmetric positive definite matrix for linear system
     * @param b inhomogeniety
     * @param if A is not symmetric positive definite throws error code 
     *
     * This function wraps the cholesky linear system solver
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
     * This function wraps the qrsolver under constraints
     */
    void linalg_constrained_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::matrix<double> &constr);

    /**
     * \brief eigenvalues of a symmetric matrix A*x=E*x
     * @param A symmetric matrix 
     * @param E vector of eigenvalues
     * @param V matrix of eigenvalues
     * 
     * This function wraps gsl_eigen_symmv / DSYEV
     * note that the eigenvalues/eigenvectors are UNSORTED 
     * 
     */
    bool linalg_eigenvalues_symmetric( ub::symmetric_matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V );
    
   /**
     * \brief eigenvalues of a symmetric matrix A*x=E*x
     * @param A matrix 
     * @param E vector of eigenvalues
     * @param V matrix of eigenvalues
     * 
     * This function wraps gsl_eigen_symmv / DSYEV
     * 
     */
    bool linalg_eigenvalues( ub::matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V );
    
    
   /**
     * \brief eigenvalues of a symmetric matrix A*x=E*x
     * @param E vector of eigenvalues
     * @param V input: matrix to diagonalize
     * @param V output: eigenvectors      
     * 
     * This function wrapps gsl_eigen_symmv / DSYEV
     * 
     */
    bool linalg_eigenvalues( ub::vector<double> &E, ub::matrix<double> &V );
    
    
       /**
     * \brief eigenvalues of a symmetric matrix A*x=E*x
     * @param E vector of eigenvalues
     * @param V input: matrix to diagonalize
     * @param V output: eigenvectors      
     * 
     * This function wrapps gsl_eigen_symmv / DSYEV
     * 
     */
    bool linalg_eigenvalues( ub::vector<float> &E, ub::matrix<float> &V );
    
   /**
     * \brief eigenvalues of a symmetric matrix A*x=E*x
     * @param E vector of eigenvalues
     * @param V input: matrix to diagonalize
     * @param V output: eigenvectors      
     * 
     * This function wrapps gsl_eigen_symmv / DSYEV
     * 
     */
    bool linalg_eigenvalues( ub::matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V , int nmax );
    
      /**
     * \brief eigenvalues of a symmetric matrix A*x=E*x single precision
     * @param E vector of eigenvalues
     * @param V input: matrix to diagonalize
     * @param V output: eigenvectors      
     * 
     * This function wrapps gsl_eigen_symmv / DSYEV
     * 
     */
    bool linalg_eigenvalues( ub::matrix<float> &A, ub::vector<float> &E, ub::matrix<float> &V , int nmax );
    
     /**
     * \brief eigenvalues of a symmetric matrix A*x=E*B*x double precision
     * @param E vector of eigenvalues
     * @param A input: matrix to diagonalize
     * @param B input: overlap matrix
     * @param V output: eigenvectors      
     * 
     * This function wrapps gsl_eigen_gensymmv / dsygv
     * 
     */
    bool linalg_eigenvalues_general( ub::matrix<double> &A,ub::matrix<double> &B, ub::vector<double> &E, ub::matrix<double> &V);
    /**
     * \brief computes Singular value decomposition A = U S V^T double precision
     * @param S vector of singular values
     * @param A input: matrix for SVD, is oerwritten with U
     * @param B input: overlap matrix
     * @param V output: eigenvectors      
     * 
     * This function wrapps eigen_gensymmv / dsygv
     * 
     */
   bool linalg_singular_value_decomposition(ub::matrix<double> &A, ub::matrix<double> &VT, ub::vector<double> &S );
    
   /**
     * \brief inverts A via svd
     * @param A symmetric positive definite matrix
     * @param V inverse matrix
     * @param lower limit of condition number of the matrix, singular values below that will be set to zero
     * This function wraps the inversion of a matrix via svd
     */
   int linalg_invert_svd(ub::matrix<double> &A, ub::matrix<double> &V,double limitCN);
   
   /**
     * \brief calculates loewdin transformation of matrices
     * @param J matrix to transform, returns transformed matrix
     * @param S, overlap matrix, returns S-1/2
     * @param returns smallest eigenvalue of S
     * This function calculates the loewdin transformation of a matrix
     */
   double linalg_loewdin(ub::matrix<double> &J, ub::matrix<double> &S);
/**
     * \brief calculates matrix sqrt of a matrix
     * @param matrix to calculate sqrt of S, return S1/2
   
     * This function calculates the sqrt of a matrix
     */
   int linalg_matrixsqrt(ub::matrix<double> &S);
   /**
     * \brief returns the the element with the largest absolute value of a matrix
     * @param matrix to find largest value of
   
     * returns the the element with the largest absolute value of a matrix
     */
   
   
   double linalg_getMax( const ub::matrix<double>& _matrix );
   /**
    *  * \brief returns the rms value of a matrix
     * @param matrix to find  rms value of
   
     * returns the rms value of a matrix
     */
   double linalg_getRMS(const ub::matrix<double>& _matrix );
   
   
   
}}



#endif	/* __VOTCA_TOOLS_LINALG_H */

