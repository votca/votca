#ifndef _GSL_BOOST_UBLAS_MATRIX_PROD_
#define _GSL_BOOST_UBLAS_MATRIX_PROD_

#ifdef NDEBUG

#include <boost/version.hpp>
#if defined (BOOST_VERSION) && (BOOST_VERSION >= 103401)

#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_blas.h>

namespace boost { namespace numeric { namespace ublas {

    //template<class T>
    matrix<double>    // prod( m1, m2 )
    prod(matrix<double> &m1, matrix<double> &m2)
    {
       gsl_matrix_view mA = gsl_matrix_view_array (&m1(0,0), m1.size1(), m1.size2());
       gsl_matrix_view mB = gsl_matrix_view_array (&m2(0,0), m2.size1(), m2.size2());

       boost::numeric::ublas::matrix<double> AxB( m1.size1(), m2.size2() );
       gsl_matrix_view mC = gsl_matrix_view_array (&AxB(0,0), AxB.size1(), AxB.size2());

       gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &mA.matrix, &mB.matrix,
                  0.0, &mC.matrix);
       return AxB;
    }

    matrix<float>    // prod( m1, m2 )
    prod(matrix<float> &m1, matrix<float> &m2)
    {
       gsl_matrix_float_view mA = gsl_matrix_float_view_array (&m1(0,0), m1.size1(), m1.size2());
       gsl_matrix_float_view mB = gsl_matrix_float_view_array (&m2(0,0), m2.size1(), m2.size2());

       boost::numeric::ublas::matrix<float> AxB( m1.size1(), m2.size2() );
       gsl_matrix_float_view mC = gsl_matrix_float_view_array (&AxB(0,0), AxB.size1(), AxB.size2());

       gsl_blas_sgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &mA.matrix, &mB.matrix,
                  0.0, &mC.matrix);
       return AxB;
    }

}}}

#endif  // BOOST_VERSION
#endif  // NDEBUG
#endif  // _GSL_BOOST_UBLAS_MATRIX_PROD_   
