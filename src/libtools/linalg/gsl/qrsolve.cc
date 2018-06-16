/* 
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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
#include <boost/numeric/ublas/io.hpp>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>


namespace votca { namespace tools {

using namespace std;


void linalg_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::vector<double> *residual)
{
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

    gsl_matrix_view m
        = gsl_matrix_view_array (&A(0,0), A.size1(), A.size2());

    gsl_vector_view gb
        = gsl_vector_view_array (&b(0), b.size());

    gsl_vector *gsl_x = gsl_vector_alloc (x.size());
    gsl_vector *tau = gsl_vector_alloc (x.size());
    gsl_vector *gsl_residual = gsl_vector_alloc (b.size());

    gsl_linalg_QR_decomp (&m.matrix, tau);

    gsl_linalg_QR_lssolve (&m.matrix, tau, &gb.vector, gsl_x, gsl_residual);

    for (size_t i =0 ; i < x.size(); i++)
        x(i) = gsl_vector_get(gsl_x, i);

    if(residual)
        for (size_t i =0 ; i < residual->size(); i++)
            (*residual)(i) = gsl_vector_get(gsl_residual, i);

    gsl_vector_free (gsl_x);
    gsl_vector_free (tau);
    gsl_vector_free (gsl_residual);
}


void linalg_constrained_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::matrix<double> &B)
{

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
            throw std::runtime_error("constrained_qrsolve_zero_column_in_matrix");
        }
    }

    const int NoEquations = b.size();
    const int NoUnknowns = x.size();
    const int Noconstrains = B.size1(); //number of constraints is number of rows of constr
    
    // Transpose constr:
    B = ub::trans(B);

    // temporary variables
    ub::matrix<double> Q(NoUnknowns, NoUnknowns);  
    // Q matrix: QR decomposition of trans(B)
    ub::matrix<double> Q_k(NoUnknowns, NoUnknowns);
    ub::identity_matrix<double> I (NoUnknowns);
    ub::vector<double> v(NoUnknowns);

    Q = ub::zero_matrix<double>(NoUnknowns, NoUnknowns);
    Q_k = ub::zero_matrix<double>(NoUnknowns, NoUnknowns);
    v = ub::zero_vector<double>(NoUnknowns);

    double *tmp = & B(0,0);
    gsl_matrix_view gsl_constr
      = gsl_matrix_view_array (tmp, B.size1(), B.size2());

    tmp = &b(0);
    gsl_vector_view gsl_b
         = gsl_vector_view_array (tmp, b.size());


    //vector to store Householder coefficients of QR-decomposition of trans(B)
    gsl_vector *tau_qr = gsl_vector_alloc (Noconstrains);
    gsl_linalg_QR_decomp (&gsl_constr.matrix, tau_qr);
    //construct the matrix Q of trans(B)=Q * R
    Q = I;

    //k=1,..,ysize (number of rows of B, meaning number of constraints)
    for (int k = Noconstrains; k > 0 ; k--) {

        for (int icout = 0; icout < k - 1; icout++) {
             v(icout) = 0;
        }
        v(k - 1) = 1.0;

        for (int icout = k; icout < NoUnknowns; icout++) {
             v(icout) = gsl_matrix_get(&gsl_constr.matrix, icout, k - 1 );
        }

        Q_k = I - gsl_vector_get(tau_qr, k - 1 ) * ub::outer_prod ( v, v );
        Q = ub::prec_prod(Q, Q_k);

    }
    //now Q is the one of trans(B)=Q * R
    Q = ub::trans(Q);
    gsl_vector_free (tau_qr);
    // Calculate A * Q and store the result in A
    A = ub::prec_prod(A, Q);

    // A = [A1 A2], so A2 is just a block of A
    // [A1 A2] has N rows. A1 has ysize columns (range(0,ysize) ).
    //A2 has 2*ngrid-ysize columns ( range(ysize,2*ngrid) ).
    ub::matrix<double> A2 = ub::matrix_range<ub::matrix<double> >(A,
            ub::range (0, NoEquations), ub::range (Noconstrains, NoUnknowns)
         );

    tmp = &A2(0,0);
    gsl_matrix_view gsl_A2
         = gsl_matrix_view_array (tmp, A2.size1(), A2.size2());
    

    //now perform QR-decomposition of A2 to solve the least-squares problem A2 * z = b
    //A2 has N rows and (2*ngrid-ysize) columns -> 
    gsl_vector *z = gsl_vector_alloc (NoUnknowns-Noconstrains);
    gsl_vector *tau_solve = gsl_vector_alloc (NoUnknowns-Noconstrains);  
    gsl_vector *residual = gsl_vector_alloc (NoEquations);
    gsl_linalg_QR_decomp (&gsl_A2.matrix, tau_solve);
    gsl_linalg_QR_lssolve (&gsl_A2.matrix, tau_solve, &gsl_b.vector, z, residual);
    // Next two cycles assemble vector from y (which is zero-vector) and z
    // (which we just got by gsl_linalg_QR_lssolve)

    //vector y has ysize components
    for (int i = 0; i < Noconstrains; i++ ) {
           x[i] = 0.0;
    }

    for (int i = Noconstrains; i < NoUnknowns; i++ ) {
           x[i] = gsl_vector_get(z, i - Noconstrains);
    }
    
    x = ub::prec_prod( Q, x );
    
    gsl_vector_free (z);
    gsl_vector_free (tau_solve);
    gsl_vector_free (residual);
}

}}
