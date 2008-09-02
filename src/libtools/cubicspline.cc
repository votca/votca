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
    
    ub::matrix<double> A(N + ngrid, 2*ngrid);
    A = ub::zero_matrix<double>(N + ngrid, 2*ngrid);
    ub::vector<double> b(N + ngrid);
    b  = ub::zero_vector<double>(N+ngrid);

    // construct the matrix to fit the points
    AddToFitMatrix(A, x, 0);
   
    ub::vector_range< ub::vector<double> > b1(b, ub::range(0,N));
    b1 = -y;
    
    // construct the smoothing condition
    // first derivatives have to be continuous:
    //AddBCToFitMatrix(A, N);
        
    
    // Solving linear equations system
        
    // libgsl has problems with solving overdetermined systems by means of 
    // QR decomposition if the system doesn't have the full rank.
    // MATLAB can handle this problem easily.
        
    double* pointer_m = &A(0,0);
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

