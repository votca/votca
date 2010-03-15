#include "linalg.h"
#include <gsl/gsl_linalg.h>

namespace votca { namespace tools {

void linalg_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b)
{
    // now A is a tridiagonal system, solve it!
    double* pointer_m = &A(0,0);
    double* pointer_b = &b(0);

    const int N = x.size();
    
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
        x(i) = gsl_vector_get(gx, i);
    }

    gsl_vector_free (gx);
    gsl_vector_free (tau);
    gsl_vector_free (residual);
}

}}
