/*
 *            Copyright 2009-2017 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#include <votca/xtp/linalg.h>

//this file is temporary till conversion to eigen in all of votca is realized


namespace votca { namespace xtp {
  
   
    bool linalg_eigenvalues(Eigen::MatrixXd&A, Eigen::VectorXd &E, Eigen::MatrixXd&V , int nmax ){
      #if defined(MKL)
    double wkopt;
    double* work;
    double abstol, vl, vu;
     
    MKL_INT lda;
    MKL_INT info;
    MKL_INT lwork;
    MKL_INT il, iu, m, ldz ;
    
    int n = A.rows();
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
    // make a pointer to the EIGEN matrix so that LAPACK understands it
    double * pA = A.data();   
    double * pV = V.data();
    double * pE = E.data();
    
    // call LAPACK via C interface
    info = LAPACKE_dsyevx( LAPACK_COL_MAJOR, 'V', 'I', 'U', n, pA , lda, vl, vu, il, iu, abstol, &m, pE, pV, n,  ifail );

   
#else
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
      V =  es.eigenvectors().block(0,0,A.rows(),nmax);
      E = es.eigenvalues().segment(0,nmax);
      int info=es.info();
#endif
       if( info != 0 ) {
        return false;
    } else {
        return true;
    }
    }
    
     
    bool linalg_eigenvalues(Eigen::MatrixXf&A, Eigen::VectorXf &E, Eigen::MatrixXf&V , int nmax ){
       #if defined(MKL)
 
    float wkopt;
    float* work;
    float abstol, vl, vu;
    MKL_INT lda;
    MKL_INT info;
    MKL_INT lwork;
    MKL_INT il, iu, m, ldz ;
    
    int n = A.rows();
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
    // make a pointer to the eigen matrix so that LAPACK understands it
    float * pA = A.data();   
    float * pV = V.data();
    float * pE = E.data();
    // call LAPACK via C interface
    info = LAPACKE_ssyevx( LAPACK_COL_MAJOR, 'V', 'I', 'U', n, pA , lda, vl, vu, il, iu, abstol, &m, pE, pV, n,  ifail );

   
#else
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(A);
      V =  es.eigenvectors().block(0,0,A.rows(),nmax);
      E = es.eigenvalues().segment(0,nmax);
      int info=es.info();
#endif
    
 if( info != 0 ) {
        return false;
    } else {
        return true;
    }
    }

}}
