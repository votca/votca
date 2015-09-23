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


void linalg_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::vector<double> *residual){
    // check matrix for zero column
    int nonzero_found = 0;
    for(size_t j=0; j<A.size2(); j++) {
        nonzero_found = 0;
        for(size_t i=0; i<A.size1(); i++) {
            if(fabs(A(i,j))>0) {
                nonzero_found = 1;
            }
        }
        if(nonzero_found==0) {
            throw "qrsolve_zero_column_in_matrix";
        }
    }

    
    MKL_INT info;
    MKL_INT sizeM = A.size1();
    MKL_INT sizeN = A.size2();
    MKL_INT nrhs = 1;
    char trans = 'N';
    
    // pointer for LAPACK
    double * pA = const_cast<double*>(&A.data().begin()[0]);
    double * pb = const_cast<double*>(&b.data()[0]);


    
    info = LAPACKE_dgels( LAPACK_ROW_MAJOR , trans , sizeM, sizeN, nrhs , pA , sizeM , pb , sizeM );

    if ( info != 0 ) 
        throw std::runtime_error("QR least-squares solver failed");


    for (size_t i =0 ; i < x.size(); i++){
        x(i) = b(i);
    }
    
    for (size_t i = 0 ; i < b.size(); i++){
        (*residual)(i) = b(i + x.size() );
    }
}



void linalg_constrained_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::matrix<double> &constr){
        // matrix inversion using MKL
    throw std::runtime_error("linalg_constrained_qrsolve is not compiled-in due to disabling of GSL - recompile Votca Tools with GSLsupport");
}

}}
