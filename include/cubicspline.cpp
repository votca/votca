/* 
 * File:   cubicspline.cpp
 * Author: ruehle
 *
 * Created on August 22, 2008, 4:44 PM
 */

#include <cubicspline.h>
#include <gsl/gsl_linalg.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>

using namespace std;

void CubicSpline::Fit(ub::vector<double> x, ub::vector<double> y)
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
    
    ub::matrix<double> A(N + ngrid, 2*ngrid);
    ub::vector<double> b(2*ngrid);
    
    // construct the matrix to fit the points
    for(int i=0; i<N; ++i) {
        int spi = getInterval(x[i]);
        A(i, spi) = this->A(x[i]);
        A(i, spi+1) = B(x[i]);
        A(i, spi + ngrid) = C(x[i]);
        A(i, spi + ngrid + 1) = D(x[i]);
        b(i) = y(i);
    }
    
    // construct the smoothing condition
    // first derivatives have to be continuous:
    // A'_i(x[i+1])*f[i] + B'_i(x[i+1])*f[i+1] + C'_i(x[i+1])*f''[i] + D'_i(x[i+1])*f''[i+1]
    // = A'_{i+1}(x[i+1])*f[i] + B'_{i+1}(x[i+2])*f[i+1] + C'_{i+1}(x[i+1])*f''[i] + D'_{i+i}(x[i+1])*f''[i+2]
    for(int i=0; i<ngrid - 2; ++i) {
            A(N+i,i) = A_prime(i, 0);
            A(N+i, i+1) = B_prime(i,0) - A_prime(i,1);
            A(N+i, i+2) = -B_prime(i,1);

            A(N+i, ngrid+1 + i) = C_prime(i,0);
            A(N+i, ngrid+1 + i+1) = D_prime(i,0) - C_prime(i,1);
            A(N+i, ngrid+1 + i+2) = -D_prime(i,1);
    }
    
    // currently only natural boundary conditions:
    A(N + ngrid - 2, ngrid) = 0;
    A(N + ngrid - 1, 2*ngrid-1) = 0;
 
    
            // Solving linear equations system
        
        // libgsl has problems with solving overdetermined systems by means of 
        // QR decomposition if the system doesn't have the full rank.
        // MATLAB can handle this problem easily.
        
    double* pointer_m = &A(0,0);
    double* pointer_b = & b(0);        

    gsl_matrix_view m
        = gsl_matrix_view_array (pointer_m, b.size(), 2*(ngrid+1));
    
    gsl_vector_view gb
        = gsl_vector_view_array (pointer_b, b.size());
    
    gsl_vector *gx = gsl_vector_alloc (2*(ngrid+1));
    gsl_vector *tau = gsl_vector_alloc (2*(ngrid+1));       
    gsl_vector *residual = gsl_vector_alloc (b.size());
    
    gsl_linalg_QR_decomp (&m.matrix, tau);

    gsl_linalg_QR_lssolve (&m.matrix, tau, &gb.vector, gx, residual);
             
    std::cout << "write out results\n";
    for (int i =0 ; i < 2*(ngrid+1); i++) {
        cout << "x[" << i <<"] = " << gsl_vector_get(gx, i) << "\n";
        _r(i) = gsl_vector_get(gx, i);
    }
        
    cout << "r  " << "F(r)  " << endl;
        
    PrintOutResult();
        
    gsl_vector_free (gx);
    gsl_vector_free (tau);
    gsl_vector_free (residual);
}

