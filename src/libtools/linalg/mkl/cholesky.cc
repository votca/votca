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

#include <votca/tools/linalg.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <votca/tools/votca_config.h>

#ifndef NOMKL
#include "mkl.h"
#include "mkl_lapacke.h"
#endif


namespace votca { namespace tools {

using namespace std;


void linalg_cholesky_decompose( ub::matrix<double> &A){
    
#ifdef NOMKL
    throw std::runtime_error("linalg_cholesky_decompose is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
#else
    // Cholesky decomposition using MKL
    // input matrix A will be changed

    // LAPACK variables
    MKL_INT info;
    MKL_INT n = A.size1();
    char uplo = 'L';
    
    // pointer for LAPACK
    double * pA = const_cast<double*>(&A.data().begin()[0]);
    info = LAPACKE_dpotrf( LAPACK_ROW_MAJOR , uplo , n, pA, n );
    if ( info != 0 )
        throw std::runtime_error("Matrix not symmetric positive definite");
    
#endif   
}




void linalg_cholesky_solve( ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b ){
    
#ifdef NOMKL
    throw std::runtime_error("linalg_cholesky_solve is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
#else
      
    /* calling program should catch the error error code
     * thrown by LAPACKE_dpotrf and take
     * necessary steps
     */
    
    
    // LAPACK variables
    MKL_INT info;
    MKL_INT n = A.size1();
    char uplo = 'L';
    
    // pointer for LAPACK LU factorization of input matrix
    double * pA = const_cast<double*>(&A.data().begin()[0]); // input array
     
    // get LU factorization
    info = LAPACKE_dpotrf( LAPACK_ROW_MAJOR , uplo , n, pA, n );
    
    if ( info != 0 )
        throw std::runtime_error("Matrix not symmetric positive definite");
    
    MKL_INT nrhs = 1;
    
    // pointer of LAPACK LU solver
    double * pb = const_cast<double*>(&b.data()[0]);
    info = LAPACKE_dpotrs(LAPACK_ROW_MAJOR, uplo, n, nrhs, pA, n, pb, n );

    // on output, b contains solution
    x = b;

    
#endif   
}




}}
