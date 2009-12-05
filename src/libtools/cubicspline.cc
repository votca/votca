/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

/* 
 * File:   cubicspline.cpp
 * Author: ruehle
 *
 * Created on August 22, 2008, 4:44 PM
 */

#include <cubicspline.h>
#include <gsl/gsl_linalg.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
//#include <householderqr.h>
#include <iostream>

using namespace std;

void CubicSpline::Interpolate(ub::vector<double> &x, ub::vector<double> &y)
{    
    if(x.size() != y.size())
        throw std::invalid_argument("error in CubicSpline::Fit : size of vector x and y dos not match");
    
    const int N = x.size();
    
    // adjust the grid
    _r.resize(N);
    _f.resize(2*N);    
    
    // create vector proxies to individually access f and f''
    ub::vector_range<ub::vector<double> > f (_f, ub::range (0, N));
    ub::vector_range<ub::vector<double> > f2 (_f, ub::range (N, 2*N));

    // copy the grid points into f
    _r = x;
    f = y;
    f2 = ub::zero_vector<double>(N);
    
    // not calculate the f''
    ub::matrix<double> A(N, N);
    A = ub::zero_matrix<double>(N,N);
    
    for(int i=0; i<N - 2; ++i) {
            f2(i+1) = -( A_prime_l(i)*f(i)
            + (B_prime_l(i) - A_prime_r(i)) * f(i+1)
            -B_prime_r(i) * f(i+2));

            A(i+1, i) = C_prime_l(i);
            A(i+1, i+1) = D_prime_l(i) - C_prime_r(i);
            A(i+1, i+2) = -D_prime_r(i);
    }
    
    switch(_boundaries) {
        case splineNormal:
            A(0, 0) = 1;
            A(N - 1, N-1) = 1;
            break;
        case splinePeriodic:
            A(0,0) = 1; A(0,N-1) = -1;
            A(N-1,0) = 1; A(N-1,N-1) = -1;
            break;
    }
    
    

    // now A is a tridiagonal system, solve it!
    double* pointer_m = &A(0,0);
    double* pointer_b = &f2(0);        

    gsl_matrix_view m
        = gsl_matrix_view_array (pointer_m, N, N);
    
    gsl_vector_view gb
        = gsl_vector_view_array (pointer_b, N);
    
    gsl_vector *gx = gsl_vector_alloc (N);
    gsl_vector *tau = gsl_vector_alloc (N);       
    gsl_vector *residual = gsl_vector_alloc (N);
    
    gsl_linalg_QR_decomp (&m.matrix, tau);

    gsl_linalg_QR_lssolve (&m.matrix, tau, &gb.vector, gx, residual);
             
    for (int i =0 ; i < N; i++) {
        f2(i) = gsl_vector_get(gx, i);
    }
        
    gsl_vector_free (gx);
    gsl_vector_free (tau);
    gsl_vector_free (residual);
}

void CubicSpline::Fit(ub::vector<double> &x, ub::vector<double> &y)
{
    if(x.size() != y.size())
        throw std::invalid_argument("error in CubicSpline::Fit : size of vector x and y dos not match");
    
    const int N = x.size();
    const int ngrid = _r.size();
    
    // construct the equation
    // A*u = b
    // where u = { {f[i]}, {f''[i]} }
    // and b[i] = y[i] for 0<=i<N
    // and b[i]=0 for i>=N (for smoothing condition)
    // A[i,j] contains the data fitting + the spline smoothing conditions
    
    ub::matrix<double> A(N, 2*ngrid);
    A = ub::zero_matrix<double>(N, 2*ngrid);
    ub::vector<double> b(N);
    b  = ub::zero_vector<double>(N);
    
    ub::matrix<double> B_constr(ngrid, 2*ngrid);  // Matrix with smoothing conditions
    B_constr = ub::zero_matrix<double>(ngrid, 2*ngrid);
    
    ub::matrix<double> Q(2*ngrid, 2*ngrid);       // Q matrix: QR decomposition of trans(B)
    Q = ub::zero_matrix<double>(2*ngrid, 2*ngrid);
    
    ub::matrix<double> A2(N, ngrid);             // Matrix A2 (see manual)
    A2 = ub::zero_matrix<double>(N, ngrid);
    
    // Construct smoothing matrix
    AddBCToFitMatrix(B_constr, 0);
    
    // Transpose B_constr:
    B_constr = trans(B_constr);
    
//================ START QR =============================//
    ub::matrix<double> Q_k(2*ngrid, 2*ngrid);
    Q_k = ub::zero_matrix<double>(2*ngrid, 2*ngrid);
       
    ub::identity_matrix<double> I (2*ngrid); 

    ub::vector<double> v(2*ngrid); 
    v = ub::zero_vector<double>(2*ngrid);
    
    double* pointer_m = & B_constr(0,0);    
    
    gsl_matrix_view B_t 
      = gsl_matrix_view_array (pointer_m, 2*ngrid, ngrid);
     
    gsl_vector *tau_qr = gsl_vector_alloc (ngrid); 
    
    gsl_linalg_QR_decomp (&B_t.matrix, tau_qr);   

    Q = I;
              
    for (int k = ngrid; k > 0 ; k--) {
           
        for (int icout = 0; icout < k - 1; icout++) {
             v(icout) = 0;
        }
        v(k - 1) = 1.0;

        for (int icout = k; icout < 2*ngrid; icout++) {
             v(icout) = gsl_matrix_get(&B_t.matrix, icout, k - 1 );
        }
           
        Q_k = I - gsl_vector_get(tau_qr, k - 1 ) * outer_prod ( v, v );
        Q = prec_prod(Q, Q_k);
           
    }
    
    Q = trans(Q);  
    gsl_vector_free (tau_qr);
    
//=============== END QR ====================================//
    
    // construct the matrix to fit the points and the vector b
    AddToFitMatrix(A, x, 0);
    b = -y; // why is it -y?
    
    // Calculate A * Q and store the result in A
    A = prec_prod(A, Q); 
    
    // A = [A1 A2], so A2 is just a block of A
    for (int iraw = 0; iraw < N; iraw++) {
        for (int icol = ngrid; icol < 2*ngrid; icol++) {
            A2(iraw, icol - ngrid) = A(iraw, icol);
        }
    }    
    
//============== START SOLVE QR-DECOMPOSITION================//
    pointer_m = & A2(0,0);
    double* pointer_b = & b(0);
    
    gsl_matrix_view m
         = gsl_matrix_view_array (pointer_m, N, ngrid);
    
    gsl_vector_view gsl_b
         = gsl_vector_view_array (pointer_b, N);
    
    gsl_vector *z = gsl_vector_alloc (ngrid);
  gsl_vector *tau_solve = gsl_vector_alloc (ngrid);  // already done!         
    gsl_vector *residual = gsl_vector_alloc (N);    
    
    gsl_linalg_QR_decomp (&m.matrix, tau_solve);
    gsl_linalg_QR_lssolve (&m.matrix, tau_solve, &gsl_b.vector, z, residual);
    
    // Next two cycles assemble vector from y (which is zero-vector) and z 
    // (which we just got by gsl_linalg_QR_lssolve)
    
    for (int i = 0; i < ngrid; i++ ) {
           _f[i] = 0.0;
    }  
    
    for (int i = ngrid; i < 2 * ngrid; i++ ) {
           _f[i] = gsl_vector_get(z, i - ngrid);
    }      
    
    // To get the final answer this vector should be multiplied by matrix Q
    // TODO: here i changed the sign, check again! (victor)
    _f = -prec_prod( Q, _f );    
    
    gsl_vector_free (z);
    gsl_vector_free (tau_solve);
    gsl_vector_free (residual);
      
//============== END SOLVE QR ===============================//
    
    // construct the smoothing condition
    // first derivatives have to be continuous:
    //AddBCToFitMatrix(A, N);
        
    
    // Solving linear equations system
        
    // libgsl has problems with solving overdetermined systems by means of 
    // QR decomposition if the system doesn't have the full rank.
    // MATLAB can handle this problem easily.
        
/*    double* pointer_m = &A(0,0);
    double* pointer_b = & b(0);        

    gsl_matrix_view m
        = gsl_matrix_view_array (pointer_m, N + ngrid, 2*ngrid);
    
    gsl_vector_view gb
        = gsl_vector_view_array (pointer_b, N+ngrid);
    
    gsl_vector *gx = gsl_vector_alloc (2*ngrid);
    gsl_vector *tau = gsl_vector_alloc (2*ngrid);       
    gsl_vector *residual = gsl_vector_alloc (b.size());
    
    gsl_linalg_QR_decomp (&m.matrix, tau);

    gsl_linalg_QR_lssolve (&m.matrix, tau, &gb.vector, gx, residual);
             
    for (int i =0 ; i < 2*ngrid; i++) {
        _f(i) = gsl_vector_get(gx, i);
    }
        
    gsl_vector_free (gx);
    gsl_vector_free (tau);
    gsl_vector_free (residual);
  */  
    return;
    
  /*  using namespace boost::numeric::ublas;
    using namespace std;

    matrix<double> Q(3,3), R(3,3);
    HouseholderQR? (A,Q,R);
    matrix<double> Z = prod(Q,R) - A;
    float f = norm_1 (Z);
    cout << "Q=" << Q <<endl;
    cout << "R=" << R << endl;
    cout << "|Q*R - A|=" << f << endl;*/

}

