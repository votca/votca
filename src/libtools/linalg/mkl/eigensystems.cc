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


/* Wrapper for DSYEV
 *  - make copy of input matrix
 */


bool linalg_eigenvalues( ub::matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V)
{
#ifdef NOMKL
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
#else
    
    // cout << " \n I'm really using MKL! " << endl;
    
    int n = A.size1();
    int lda = n ;
    // make sure that containers for eigenvalues and eigenvectors are of correct size
    E.resize(n);
    V.resize(n, n);
    // Query and allocate the optimal workspace 
    double wkopt;
    double* work;
    int info;
    int lwork;
    lwork = -1;

    // MKL is different to GSL because it overwrites the input matrix
    V = A; // make a copy (might actually be unnecessary in most cases!)

    // make a pointer to the ublas matrix so that LAPACK understands it
    double * pV = const_cast<double*>(&V.data().begin()[0]);
    double * pE = const_cast<double*>(&E.data()[0]);
    
    // call LAPACK via C interface
    info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, pV , lda, pE );

    if( info > 0 ) {
        return false;
    } else {
        return true;
    }

#endif
};



bool linalg_eigenvalues(  ub::vector<double> &E, ub::matrix<double> &V)
{
#ifdef NOMKL
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
#else
    
    // cout << " \n I'm really using MKL! " << endl;
    
    int n = V.size1();
    int lda = n ;
    // make sure that containers for eigenvalues and eigenvectors are of correct size
    E.resize(n);
    // V.resize(n, n);
    // Query and allocate the optimal workspace 
    double wkopt;
    double* work;
    int info;
    int lwork;
    lwork = -1;

    // MKL is different to GSL because it overwrites the input matrix
    // V = A; // make a copy (might actually be unnecessary in most cases!)

    // make a pointer to the ublas matrix so that LAPACK understands it
    double * pV = const_cast<double*>(&V.data().begin()[0]);
    double * pE = const_cast<double*>(&E.data()[0]);
    
    // call LAPACK via C interface
    info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, pV , lda, pE );

    if( info > 0 ) {
        return false;
    } else {
        return true;
    }

#endif
};





}}
