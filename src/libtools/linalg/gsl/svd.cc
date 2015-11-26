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
/*
bool linalg_singular_value_decomposition(ub::matrix<double> &A, ub::matrix<double> &VT, ub::vector<double> &S )
{
	
        gsl_error_handler_t *handler = gsl_set_error_handler_off();
	const size_t M = A.size1();
        const size_t N = A.size2();
        // gsl does not handle conversion of a symmetric_matrix 
        
        
	S.resize(N, false);
	VT.resize(M, N, false);
        
	gsl_matrix_view A_view = gsl_matrix_view_array(&A(0,0), M, N);
	gsl_vector_view S_view = gsl_vector_view_array(&S(0), N);
	gsl_matrix_view V_view = gsl_matrix_view_array(&V(0,0), M, N);
	gsl_vector * work = gsl_vector_alloc(N);

        int status = gsl_linalg_SV_decomp (gsl_matrix * A, gsl_matrix * V, gsl_vector * S, gsl_vector * work)
	//gsl_eigen_symmv_sort(&E_view.vector, &V_view.matrix, GSL_EIGEN_SORT_ABS_ASC);
	gsl_set_error_handler(handler);
        gsl_vector_free (work);
	return (status != 0);
         
    
};
*/

}}
