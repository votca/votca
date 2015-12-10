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


bool linalg_singular_value_decomposition(ub::matrix<double> &A, ub::matrix<double> &VT, ub::vector<double> &S ){
        // matrix inversion using MKL
    
    
    
    // define LAPACK variables
    MKL_INT m = A.size1();
    MKL_INT n = A.size2();
    //MKL_INT info;
    //MKL_INT ipiv[n];
    ub::matrix<double>work=ub::zero_matrix<double>(m,n);
    // initialize V
    S.resize(n, false);
    VT.resize(n, n, false);
    
    // pointers for LAPACK
    double * a = const_cast<double*>(&A.data().begin()[0]);
    double * s = const_cast<double*>(&S.data().begin()[0]);
    double * vt = const_cast<double*>(&VT.data().begin()[0]);   
    double * superb = const_cast<double*>(&work.data().begin()[0]);
    // solve
    int status= LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'O', 'A',  m,  n,  a, n,  s,  NULL,m,  vt, n, superb );
    
    return (status != 0);
    
}

}}
