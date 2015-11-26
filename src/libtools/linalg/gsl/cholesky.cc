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

void linalg_cholesky_decompose( ub::matrix<double> &A){
        // Cholesky decomposition using GSL
        const size_t N = A.size1();
        
        gsl_matrix_view A_view = gsl_matrix_view_array(&A(0,0), N, N);
        
        // get the Cholesky matrices
        (void)gsl_linalg_cholesky_decomp ( &A_view.matrix );
}

void linalg_cholesky_decompose( ub::matrix<float> &A){
        // Cholesky decomposition using GSL
    throw std::runtime_error("linalg_cholesky_decompose (float) is not compiled-in due to disabling of MKL - recompile Votca Tools with MKL support");
   
}




void linalg_cholesky_solve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b){
    /* calling program should catch the error error code GSL_EDOM
     * thrown by gsl_linalg_cholesky_decomp and take
     * necessary steps
     */
    
    gsl_matrix_view m
        = gsl_matrix_view_array (&A(0,0), A.size1(), A.size2());

    gsl_vector_view gb
        = gsl_vector_view_array (&b(0), b.size());

    gsl_vector *gsl_x = gsl_vector_alloc (x.size());

    gsl_set_error_handler_off();
    int status = gsl_linalg_cholesky_decomp(&m.matrix);

    if( status == GSL_EDOM)
        throw std::runtime_error("Matrix not symmetric positive definite");

    
    gsl_linalg_cholesky_solve(&m.matrix, &gb.vector, gsl_x);

    for (size_t i =0 ; i < x.size(); i++)
        x(i) = gsl_vector_get(gsl_x, i);

    gsl_vector_free (gsl_x);
}


}}
