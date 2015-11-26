/* 
 * Copyright 2009-2015 The VOTCA Development Team (http://www.votca.org)
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

#include "mkl.h"
#include "mkl_lapacke.h"

namespace votca { namespace tools {

using namespace std;


void linalg_invert( ub::matrix<double> &A, ub::matrix<double> &V){
    // matrix inversion using MKL
    // input matrix is destroyed, make local copy
    ub::matrix<double> work = A;
    
    // define LAPACK variables
    MKL_INT n = A.size1();
    MKL_INT info;
    MKL_INT ipiv[n];
    
    // initialize V
    V = ub::identity_matrix<double>(n,n);

    // pointers for LAPACK
    double * pV = const_cast<double*>(&V.data().begin()[0]);
    double * pwork = const_cast<double*>(&work.data().begin()[0]);
    
    // solve
    info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, n, pwork , n, ipiv, pV, n );
}


void linalg_invert( ub::matrix<float> &A, ub::matrix<float> &V){
    // matrix inversion using MKL
    // input matrix is destroyed, make local copy
    ub::matrix<float> work = A;
    
    // define LAPACK variables
    MKL_INT n = A.size1();
    MKL_INT info;
    MKL_INT ipiv[n];
    
    // initialize V
    V = ub::identity_matrix<float>(n,n);

    // pointers for LAPACK
    float * pV = const_cast<float*>(&V.data().begin()[0]);
    float * pwork = const_cast<float*>(&work.data().begin()[0]);
    
    // solve
    info = LAPACKE_sgesv( LAPACK_ROW_MAJOR, n, n, pwork , n, ipiv, pV, n );
}



}}
