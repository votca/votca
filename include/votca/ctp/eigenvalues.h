#include <memory.h>
#include <gsl/gsl_eigen.h>
#include <boost/numeric/ublas/symmetric.hpp>

#ifndef __VOTCA_CTP_EIGENVALUES_H
#define __VOTCA_CTP_EIGENVALUES_H


using namespace std;

namespace votca { namespace tools {
    namespace ub = boost::numeric::ublas;
    using namespace std;
    
/**
*
* ublas binding for gsl_eigen_symmv
* note that the eigenvalues/eigenvectors are UNSORTED 
* 
*/
bool EigenvaluesSymmetric( ub::symmetric_matrix<double> &A, ub::vector<double> &E, ub::matrix<double> &V
)
{
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
};


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
bool SingularValueDecomposition( ub::matrix<double> &A, ub::matrix<double> &V, ub::vector<double> &S )
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
};


}}

#endif  /* __VOTCA_CTP_EIGENVALUES_H */