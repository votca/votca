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

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>


namespace votca { namespace tools {

using namespace std;




/**
 * ublas binding to GSL  Singular Value Decomposition
 * 
 * A = U S V^T
 * 
 * @param A MxN matrix do decompose. Becomes an MxN orthogonal matrix U
 * @param V NxN orthogonal square matrix
 * @param E NxN diagonal matrix of singular values
 * @return succeeded or not 
 */
bool linalg_singular_value_decomposition( ub::matrix<double> &A, ub::matrix<double> &V, ub::vector<double> &S )
{
	/*
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
         */
    throw std::runtime_error("linalg_singular_value_decomposition is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
};


}}
