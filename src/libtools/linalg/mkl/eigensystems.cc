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


/* Wrapper for DSYEV
 *  - make copy of input matrix
 */


bool linalg_eigenvalues( const ub::matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V)
{
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
};




bool linalg_eigenvalues_symmetric( const ub::symmetric_matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V)
{
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

    // MKL does not handle conversion of a symmetric_matrix 
    V = A;
    
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
};


bool linalg_eigenvalues(  ub::vector<double> &E, ub::matrix<double> &V)
{
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
};



bool linalg_eigenvalues(  ub::vector<float> &E, ub::matrix<float> &V)
{
    // cout << " \n I'm really using MKL! " << endl;
    
    int n = V.size1();
    int lda = n ;
    // make sure that containers for eigenvalues and eigenvectors are of correct size
    E.resize(n);
    // V.resize(n, n);
    // Query and allocate the optimal workspace 
    float wkopt;
    float* work;
    int info;
    int lwork;
    lwork = -1;

    // MKL is different to GSL because it overwrites the input matrix
    // V = A; // make a copy (might actually be unnecessary in most cases!)

    // make a pointer to the ublas matrix so that LAPACK understands it
    float * pV = const_cast<float*>(&V.data().begin()[0]);
    float * pE = const_cast<float*>(&E.data()[0]);
    
    // call LAPACK via C interface
    info = LAPACKE_ssyev( LAPACK_ROW_MAJOR, 'V', 'U', n, pV , lda, pE );

    if( info > 0 ) {
        return false;
    } else {
        return true;
    }
};




/*
 * use expert routine to calculate only a subrange of eigenvalues
 */
bool linalg_eigenvalues(const ub::matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V , int nmax)
{
    /*
     * INPUT:  matrix A (N,N)
     * OUTPUT: matrix V (N,NMAX)
     *         vector E (NMAX)
     */
    double wkopt;
    double* work;
    double abstol, vl, vu;
     
    MKL_INT lda;
    MKL_INT info;
    MKL_INT lwork;
    MKL_INT il, iu, m, ldz ;
    
    int n = A.size1();
    MKL_INT ifail[n];
    lda = n;
    ldz = nmax;
    
    
    
    // make sure that containers for eigenvalues and eigenvectors are of correct size
    E.resize(nmax);
    V.resize(n,nmax);

    
    lwork = -1;
    il = 1;
    iu = nmax;
    abstol = 0.0; // use default
    vl = 0.0;
    vu = 0.0;
    // make a pointer to the ublas matrix so that LAPACK understands it
    double * pA = const_cast<double*>(&A.data().begin()[0]);   
    double * pV = const_cast<double*>(&V.data().begin()[0]);
    double * pE = const_cast<double*>(&E.data()[0]);
    
    // call LAPACK via C interface
    info = LAPACKE_dsyevx( LAPACK_ROW_MAJOR, 'V', 'I', 'U', n, pA , lda, vl, vu, il, iu, abstol, &m, pE, pV, nmax,  ifail );

    if( info != 0 ) {
        return false;
    } else {
        return true;
    }
};



/*
 * use expert routine to calculate only a subrange of eigenvalues
 */
bool linalg_eigenvalues(const ub::matrix<float> &A, ub::vector<float> &E, ub::matrix<float> &V , int nmax)
{
    /*
     * INPUT:  matrix A (N,N)
     * OUTPUT: matrix V (N,NMAX)
     *         vector E (NMAX)
     */
    float wkopt;
    float* work;
    float abstol, vl, vu;
     
    MKL_INT lda;
    MKL_INT info;
    MKL_INT lwork;
    MKL_INT il, iu, m, ldz ;
    
    int n = A.size1();
    MKL_INT ifail[n];
    lda = n;
    ldz = nmax;
    
    
    
    // make sure that containers for eigenvalues and eigenvectors are of correct size
    E.resize(nmax);
    V.resize(n,nmax);

    
    lwork = -1;
    il = 1;
    iu = nmax;
    abstol = 0.0; // use default
    vl = 0.0;
    vu = 0.0;
    // make a pointer to the ublas matrix so that LAPACK understands it
    float * pA = const_cast<float*>(&A.data().begin()[0]);   
    float * pV = const_cast<float*>(&V.data().begin()[0]);
    float * pE = const_cast<float*>(&E.data()[0]);
    
    // call LAPACK via C interface
    info = LAPACKE_ssyevx( LAPACK_ROW_MAJOR, 'V', 'I', 'U', n, pA , lda, vl, vu, il, iu, abstol, &m, pE, pV, nmax,  ifail );

    if( info != 0 ) {
        return false;
    } else {
        return true;
    }
};



/* calculate the eigenvalues and vectors of the generalized eigenvalue problem */
bool linalg_eigenvalues_general(const ub::matrix<double> &A,const ub::matrix<double> &B, ub::vector<double> &E, ub::matrix<double> &V)
{
    // cout << " \n I'm really using MKL! " << endl;
    //check to see if matrices have same size
    int lda = A.size1();
    MKL_INT n = lda ;
    int ldb =B.size1();
    ub::matrix<double> _B(ldb,ldb);
    _B=B;
    if (lda!=ldb){
        cout << "Matrices A and B have not the same size"<< endl;
        exit(1);
    }
    // make sure that containers for eigenvalues and eigenvectors are of correct size
    E.resize(lda);
    V.resize(lda, lda);
    // Query and allocate the optimal workspace 
   
    MKL_INT info=0;
    MKL_INT LDA=lda;
    MKL_INT LDB=ldb;
    // MKL is different to GSL because it overwrites the input matrix
    V = A; // make a copy (might actually be unnecessary in most cases!)

    // make a pointer to the ublas matrix so that LAPACK understands it 
    double * pB = const_cast<double*>(&_B.data().begin()[0]);  
    double * pV = const_cast<double*>(&V.data().begin()[0]);
    double * pE = const_cast<double*>(&E.data()[0]);
    
    // call LAPACK via C interface
    info = LAPACKE_dsygv( LAPACK_ROW_MAJOR,1,'V', 'U', n, pV , lda,pB,ldb, pE );

    if( info != 0 ) {
        return false;
    } else {
        return true;
    }
};


}}
