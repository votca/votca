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

void linalg_constrained_qrsolve(ub::vector<double> &x, ub::matrix<double> &A, ub::vector<double> &b, ub::matrix<double> &constr)
{
    // Transpose constr:
    constr = trans(constr);

    const int N = b.size();
    const int ngrid = x.size()/2;

    // temporary variables
    ub::matrix<double> Q(2*ngrid, 2*ngrid);       // Q matrix: QR decomposition of trans(B)
    ub::matrix<double> A2(N, ngrid);              // Matrix A2 (see manual)
    ub::matrix<double> Q_k(2*ngrid, 2*ngrid);
    ub::identity_matrix<double> I (2*ngrid);
    ub::vector<double> v(2*ngrid);

    Q = ub::zero_matrix<double>(2*ngrid, 2*ngrid);
    A2 = ub::zero_matrix<double>(N, ngrid);
    Q_k = ub::zero_matrix<double>(2*ngrid, 2*ngrid);
    v = ub::zero_vector<double>(2*ngrid);

    double* pointer_m = & constr(0,0);

    gsl_matrix_view B_t
      = gsl_matrix_view_array (pointer_m, constr.size1(), constr.size2());

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

    // Calculate A * Q and store the result in A
    A = prec_prod(A, Q);

    // A = [A1 A2], so A2 is just a block of A
    for (int iraw = 0; iraw < N; iraw++) {
        for (int icol = ngrid; icol < 2*ngrid; icol++) {
            A2(iraw, icol - ngrid) = A(iraw, icol);
        }
    }

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
           x[i] = 0.0;
    }

    for (int i = ngrid; i < 2 * ngrid; i++ ) {
           x[i] = gsl_vector_get(z, i - ngrid);
    }

    // To get the final answer this vector should be multiplied by matrix Q
    // TODO: here i changed the sign, check again! (victor)
    x = -prec_prod( Q, x );

    gsl_vector_free (z);
    gsl_vector_free (tau_solve);
    gsl_vector_free (residual);
}


}}
