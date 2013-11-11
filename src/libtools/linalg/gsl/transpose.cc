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

#ifndef NOGSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#endif


namespace votca { namespace tools {

using namespace std;


void linalg_transpose( ub::matrix<double> &A, ub::matrix<double> &V){
    
#ifdef NOGSL
    throw std::runtime_error("linalg_invert is not compiled-in due to disabling of GSL - recompile Votca Tools with GSL support");
#else
        // matrix transposition using gsl
        gsl_error_handler_t *handler = gsl_set_error_handler_off();
	const size_t N = A.size1();
	const size_t M = A.size2();
        V.resize(M, N, false);
	// Define all the used matrices
        gsl_matrix_view A_view = gsl_matrix_view_array(&A(0,0), N, M);
        gsl_matrix_view V_view = gsl_matrix_view_array(&V(0,0), M, N);

        int status = gsl_matrix_transpose_memcpy (&V_view.matrix, &A_view.matrix);

        gsl_set_error_handler(handler);
        
	// return (status != 0);
    
#endif   
}

void linalg_transpose( ub::matrix<double> &A){
    
#ifdef NOGSL
    throw std::runtime_error("linalg_transpose is not compiled-in due to disabling of GSL - recompile Votca Tools with GSL support");
#else
        // matrix transposition using gsl
        // after inversion, original matrix is overwritten!
        gsl_error_handler_t *handler = gsl_set_error_handler_off();
	const size_t N = A.size1();
	// Define all the used matrices
        gsl_matrix_view A_view = gsl_matrix_view_array(&A(0,0), N, N);

        // this only works on square matrices

        int status = gsl_matrix_transpose(&A_view.matrix);

        gsl_set_error_handler(handler);
        
	// return (status != 0);
    
#endif   
}


}}
