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
#include <boost/numeric/ublas/matrix.hpp>
#include <math.h>       /* sqrt */
#ifndef NOGSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#endif


namespace votca { namespace tools {

using namespace std;
 
/**
*
* ublas binding for gsl_eigen_symmv
* note that the eigenvalues/eigenvectors are UNSORTED 
* 
*/
bool linalg_eigenvalues_symmetric( ub::symmetric_matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V)
{
#ifdef NOGSL
    throw std::runtime_error("linalg_eigenvalues_symmetric is not compiled-in due to disabling of GSL - recompile Votca Tools with GSL support");
#else
    
	gsl_error_handler_t *handler = gsl_set_error_handler_off();
	const size_t N = A.size1();
        
        // gsl does not handle conversion of a symmetric_matrix 
        ub::matrix<double> _A( N,N );
        _A = A;
        
	E.resize(N, false);
	V.resize(N, N, false);
	gsl_matrix_view A_view = gsl_matrix_view_array(&_A(0,0), N, N);
	gsl_vector_view E_view = gsl_vector_view_array(&E(0), N);
	gsl_matrix_view V_view = gsl_matrix_view_array(&V(0,0), N, N);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N);

	int status = gsl_eigen_symmv(&A_view.matrix, &E_view.vector, &V_view.matrix, w);
	//gsl_eigen_symmv_sort(&E_view.vector, &V_view.matrix, GSL_EIGEN_SORT_ABS_ASC);
	gsl_eigen_symmv_free(w);
	gsl_set_error_handler(handler);
        
	return (status != 0);
#endif
};


/**
*
* ublas binding for gsl_eigen_symmv
* input matrix type general matrix!
* wrapping gsl_eigen_symmv 
* 
*/
bool linalg_eigenvalues( ub::matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V)
{
#ifdef NOGSL
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL - recompile Votca Tools with GSL support");
#else
    
	gsl_error_handler_t *handler = gsl_set_error_handler_off();
	const size_t N = A.size1();
        
        // gsl does not handle conversion of a symmetric_matrix 
        ub::matrix<double> _A( N,N );
        _A = A;
        
	E.resize(N, false);
	V.resize(N, N, false);
	gsl_matrix_view A_view = gsl_matrix_view_array(&_A(0,0), N, N);
	gsl_vector_view E_view = gsl_vector_view_array(&E(0), N);
	gsl_matrix_view V_view = gsl_matrix_view_array(&V(0,0), N, N);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N);

	int status = gsl_eigen_symmv(&A_view.matrix, &E_view.vector, &V_view.matrix, w);
	gsl_eigen_symmv_sort(&E_view.vector, &V_view.matrix, GSL_EIGEN_SORT_VAL_ASC);
	gsl_eigen_symmv_free(w);
	gsl_set_error_handler(handler);
        
	return (status != 0);
#endif
};

/**
*
* ublas binding for gsl_eigen_symmv
* input matrix type general matrix single precision!
* wrapping gsl_eigen_symmv 
* 
*/
bool linalg_eigenvalues( ub::matrix<float> &A, ub::vector<float> &E, ub::matrix<float> &V)
{
#ifdef NOGSL
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL - recompile Votca Tools with GSL support");
#else
    //throw std::runtime_error("linalg_eigenvalues is not implemented with float precision, switch to MKL");
    //return -1;
    
	gsl_error_handler_t *handler = gsl_set_error_handler_off();
	const size_t N = A.size1();
        
        // gsl does not handle symmetric_matrix and floats, so this is super stupid
        ub::matrix<double> _A( N,N );
        _A = A;
        ub::vector<double> _E(N);
        ub::matrix<double> _V(N,N);
	gsl_matrix_view A_view = gsl_matrix_view_array(&_A(0,0), N, N);
	gsl_vector_view E_view = gsl_vector_view_array(&_E(0), N);
	gsl_matrix_view V_view = gsl_matrix_view_array(&_V(0,0), N, N);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(N);

	int status = gsl_eigen_symmv(&A_view.matrix, &E_view.vector, &V_view.matrix, w);
	gsl_eigen_symmv_sort(&E_view.vector, &V_view.matrix, GSL_EIGEN_SORT_VAL_ASC);
	gsl_eigen_symmv_free(w);
	gsl_set_error_handler(handler);

	//E.resize(N, false);
	//V.resize(N, N, false);
        E = _E;
        V = _V;
        
	return (status != 0);
     
#endif
};

bool linalg_eigenvalues(  ub::vector<float> &E, ub::matrix<float> &V)
{
#ifdef NOGSL
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL - recompile Votca Tools with GSL support");
#else
    //throw std::runtime_error("linalg_eigenvalues is not implemented with float precision, switch to MKL");
    //return -1;
    
        /* on input V is the matrix that shall be diagonalized
         * GSL does not provide an in-place routine, so we wrap 
         * gsl_eigen_symmv for compatibility
         */
    
         // make a copy of E
         ub::matrix<float> A = V;
    
         // now call wrapper for gsl_eigen_symmv
         bool status = linalg_eigenvalues( A , E, V );
	return (status != 0);
     
#endif
};


/**
*
* ublas binding for gsl_eigen_symm
* input matrix type general matrix!
* wrapping gsl_eigen_symm leaves input matrix 
* 
*/
bool linalg_eigenvalues( ub::vector<double> &E, ub::matrix<double> &V)
{
#ifdef NOGSL
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL - recompile Votca Tools with GSL support");
#else
    
        /* on input V is the matrix that shall be diagonalized
         * GSL does not provide an in-place routine, so we wrap 
         * gsl_eigen_symmv for compatibility
         */
    
         // make a copy of E
         ub::matrix<double> A = V;
    
         // now call wrapper for gsl_eigen_symmv
         bool status = linalg_eigenvalues( A , E, V );

         return status;
         
         
#endif
};


/*
 * use expert routine to calculate only a subrange of eigenvalues
 */
bool linalg_eigenvalues( ub::matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V , int nmax)
{
#ifdef NOGSL
    throw std::runtime_error("Available only if intell compiler is used and MKL installed");
#else    
    throw std::runtime_error("Available only if intell compiler is used and MKL installed");
#endif
}

/*
 * use expert routine to calculate only a subrange of eigenvalues single precision
 */
bool linalg_eigenvalues( ub::matrix<float> &A, ub::vector<float> &E, ub::matrix<float> &V , int nmax)
{
#ifdef NOGSL
    throw std::runtime_error("Available only if intell compiler is used and MKL installed");
#else    
    // throw std::runtime_error("Available only if intell compiler is used and MKL installed");

    // now call wrapper for gsl_eigen_symmv
    bool status = linalg_eigenvalues( A , E, V );

    return status;
         


#endif
}

bool linalg_eigenvalues_general( ub::matrix<double> &A,ub::matrix<double> &B, ub::vector<double> &E, ub::matrix<double> &V)
{
#ifdef NOGSL
    throw std::runtime_error("linalg_eigenvalues is not compiled-in due to disabling of GSL - recompile Votca Tools with GSL support");
#else
    
	gsl_error_handler_t *handler = gsl_set_error_handler_off();
	const size_t N = A.size1();
        
        // gsl destroys A so use copy
        ub::matrix<double> _A( N,N );
        _A = A;
        
        ub::matrix<double> _B( N,N );
        _B=B;
        
	E.resize(N, false);
	V.resize(N, N, false);
	gsl_matrix_view A_view = gsl_matrix_view_array(&_A(0,0), N, N);
        gsl_matrix_view B_view = gsl_matrix_view_array(&_B(0,0), N, N);
	gsl_vector_view E_view = gsl_vector_view_array(&E(0), N);
	gsl_matrix_view V_view = gsl_matrix_view_array(&V(0,0), N, N);
	gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(N);

	int status = gsl_eigen_gensymmv(&A_view.matrix,&B_view.matrix, &E_view.vector, &V_view.matrix, w);
	gsl_eigen_gensymmv_sort(&E_view.vector, &V_view.matrix, GSL_EIGEN_SORT_VAL_ASC);
	gsl_eigen_gensymmv_free(w);
	gsl_set_error_handler(handler);
        
      

    
        ub::matrix<double> _temp= ub::prod(B,V);
        ub::matrix<double> n=ub::prod(ub::trans(V),_temp);
      /*  
        for (int i=0;i<n.size1();i++){
          
                for (int j=0;j<n.size2();j++){
                cout <<"n("<< i << ":"<< j <<")= " <<n(i,j)<< endl;      
                }}
        
       */ 
        
        for (int i=0;i<n.size1();i++){
            ub::matrix_range<ub::matrix<double> > column=ub::subrange( V, 0, V.size2(),i, i+1 );
            cout <<"n("<< i << ":"<< i <<")= " <<n(i,i) <<":" <<sqrt(n(i,i))<< endl; 
            for (int j=0;j<column.size1();j++){
                
           
            cout <<"V("<<i<<":"<<j<<")="<<column(j,0)<< endl;
                    
              }
            column=column/sqrt(n(i,i));
  
        }
    
	return (status != 0);
#endif
};



}}
