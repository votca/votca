#ifndef _GSL_BOOST_UBLAS_MATRIX_PROD_
#define _GSL_BOOST_UBLAS_MATRIX_PROD_

#ifdef NDEBUG

#include <boost/version.hpp>
#if defined (BOOST_VERSION) && (BOOST_VERSION >= 103401)

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#define GSLDEBUG

#include <gsl/gsl_blas.h>

namespace boost { namespace numeric { namespace ublas {

    /* 
     * GSL has separate implementations for floats and doubles
     * hence we first have specializations for these two types
     * TSeparate templates for diagonal and symmetric matreces
     * and their combinations. To check  whether a specific
     * routine is called, define GSLDEBUG
     */
    
    // partial specialization for double precision (dgemm))
    template<class F, class A>
    inline matrix<double,F,A>    // prod( m1, m2 )
    prod(const matrix<double,F,A> &m1, const matrix<double,F,A> &m2)
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;34mGSL [CMD,CM]\x1b[0;39m\n" << std::flush;
       #endif
          
       gsl_matrix_const_view mA = gsl_matrix_const_view_array (&m1(0,0), m1.size1(), m1.size2());
       gsl_matrix_const_view mB = gsl_matrix_const_view_array (&m2(0,0), m2.size1(), m2.size2());

       boost::numeric::ublas::matrix<double,F,A> AxB( m1.size1(), m2.size2() );
       gsl_matrix_view mC = gsl_matrix_view_array (&AxB(0,0), AxB.size1(), AxB.size2());

       gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &mA.matrix, &mB.matrix,
                  0.0, &mC.matrix);
       return AxB;
    }    
    
    // partial specialization for single precision (sgemm))
    template<class F, class A>
    inline matrix<float,F,A>    // prod( m1, m2 )
    prod(const matrix<float,F,A> &m1, const matrix<float,F,A> &m2)
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;34mGSL [CMS,CM]\x1b[0;39m\n" << std::flush;
       #endif
          
       gsl_matrix_float_const_view mA = gsl_matrix_float_const_view_array (&m1(0,0), m1.size1(), m1.size2());
       gsl_matrix_float_const_view mB = gsl_matrix_float_const_view_array (&m2(0,0), m2.size1(), m2.size2());

       boost::numeric::ublas::matrix<float,F,A> AxB( m1.size1(), m2.size2() );
       gsl_matrix_float_view mC = gsl_matrix_float_view_array (&AxB(0,0), AxB.size1(), AxB.size2());

       gsl_blas_sgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &mA.matrix, &mB.matrix,
                  0.0, &mC.matrix);
       return AxB;
    }     
    
    // full specialization   
    template<class T, class F, class A>
    inline matrix<T,F,A>    // prod( m1, m2 )
    prod(const matrix<T,F,A> &m1, const matrix<T,F,A> &m2)
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;34mGSL [CM,CM]\x1b[0;39m -> " << std::flush;
       #endif
          
       return prod(m1,m2);
    }    
    
    /// transpose products
    template<class T, class F, class A>
    inline matrix<T,F,A>    // prod( trans(m1), m2 )
    prod(const matrix_unary2<matrix<T,F,A>,scalar_identity<T> > &m1, const matrix<T,F,A> &m2)
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;32mGSL [CTM,CM]\x1b[0;39m -> " << std::flush;
       #endif
       boost::numeric::ublas::matrix<T,F,A> _m1 = m1;
       return prod(_m1,m2);
    }    
    
    template<class T, class F, class A>
    inline matrix<T,F,A>    // prod( m1, trans(m2) )
    prod(const matrix<T,F,A> &m1, const matrix_unary2<matrix<T,F,A>,scalar_identity<T> > &m2)
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;32mGSL [CM,CTM]\x1b[0;39m -> " << std::flush;
       #endif
          
       const boost::numeric::ublas::matrix<T,F,A> _m2 = m2;
       return prod(m1,_m2);
    }  
    
    template<class T, class F, class A>
    inline matrix<T,F,A>    // prod( trans(m1), trans(m2) )
    prod(const matrix_unary2<matrix<T,F,A>,scalar_identity<T> > &m1, const matrix_unary2<matrix<T,F,A>,scalar_identity<T> > &m2)    
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;32mGSL [CTM,CTM]\x1b[0;39m -> " << std::flush;
       #endif
          
       boost::numeric::ublas::matrix<T,F,A> _m1 = m1;
       boost::numeric::ublas::matrix<T,F,A> _m2 = m2;
       
       return prod(_m1,_m2);
    }  

    /// diagonal matrix
    template<class T, class F, class L, class A>
    inline matrix<T,F,A>    // prod( diagonal m1, m2 )
    prod(const diagonal_matrix<T,L,A> &m1, const matrix<T,F,A> &m2 )       
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;33mGSL [CDM, CM]\x1b[0;39m -> " << std::flush;
       #endif
       const matrix<T,F,A> _m1 = m1; 
       return prod(_m1,m2);
    }

    template<class T, class F, class L, class A>
    inline matrix<T,F,A>    // prod( m1, diagonal m2 )
    prod(const matrix<T,F,A> &m1, const diagonal_matrix<T,L,A> &m2 )       
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;33mGSL [CM, CDM]\x1b[0;39m -> " << std::flush;
       #endif
       const matrix<T,F,A> _m2 = m2; 
       return prod(m1,_m2);
    }

    template<class T, class F, class L, class A>
    inline matrix<T,F,A>    // prod( diagonal m1, diagonal m2 )
    prod(const diagonal_matrix<T,L,A> &m1, const diagonal_matrix<T,L,A> &m2 )       
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;33mGSL [CDM, CM]\x1b[0;39m -> " << std::flush;
       #endif
       const matrix<T,F,A> _m1 = m1; 
       const matrix<T,F,A> _m2 = m2; 
       return prod(_m1,_m2);
    }    
    
    template<class T, class F, class L, class A>
    inline matrix<T,F,A>    // prod( diagonal m1, transpose m2 )
    prod(const diagonal_matrix<T,L,A> &m1, const matrix_unary2<matrix<T,F,A>,scalar_identity<T> > &m2)   
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;33mGSL [CDM, CTM]\x1b[0;39m -> " << std::flush;
       #endif

       const matrix<T,F,A> _m1 = m1; 
       const matrix<T,F,A> _m2 = m2; 
       
       return prod(_m1,_m2);
       
    }

   template<class T, class F, class L, class A>
    inline matrix<T,F,A>    // prod( transpose m1, diagonal m2 )
    prod(const matrix_unary2<matrix<T,F,A>,scalar_identity<T> > &m1, const diagonal_matrix<T,L,A> &m2)   
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;33mGSL [CTM, CDM]\x1b[0;39m -> " << std::flush;
       #endif

       const matrix<T,F,A> _m1 = m1; 
       const matrix<T,F,A> _m2 = m2; 
       
       return prod(_m1,_m2);
       
    }
    
    // symmetric matrix 
    template<class T, class F, class A, class TRI, class L>
    inline matrix<T,F,A>    // prod( symmetric m1, m2 )
    prod(symmetric_matrix<T, TRI, L, A> &m1, matrix<T,F,A> &m2 )   
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;33mGSL [CSM, CM]\x1b[0;39m -> " << std::flush;
       #endif
          
       assert( m1.size1() == m2.size1() );
       const matrix<T,F,A> _m1 = m1; 
       
       return prod(_m1, m2 );
    }   

    template<class T, class F, class A, class TRI, class L>
    inline matrix<T,F,A>    // prod( m1, symmetric m2 )
    prod( matrix<T,F,A> &m1, symmetric_matrix<T, TRI, L, A> &m2 )   
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;33mGSL [CSM, CM]\x1b[0;39m -> " << std::flush;
       #endif
          
       assert( m1.size1() == m2.size1() );
       const matrix<T,F,A> _m2 = m2; 
       
       return prod(m1, _m2 );
    } 

    template<class T, class F, class A, class TRI, class L>
    inline matrix<T,F,A>    // prod( symmetric m1, symmetric m2 )
    prod(symmetric_matrix<T, TRI, L, A> &m1, symmetric_matrix<T, TRI, L, A> &m2 )   
    {
       #ifdef GSLDEBUG 
          std::cout << "\x1b[0;33mGSL [CSM, CM]\x1b[0;39m -> " << std::flush;
       #endif
          
       assert( m1.size1() == m2.size1() );
       const matrix<T,F,A> _m1 = m1; 
       const matrix<T,F,A> _m2 = m2; 
       
       return prod(_m1, _m2 );
    } 
        
}}}

#endif  // BOOST_VERSION
#endif  // NDEBUG
#endif  // _GSL_BOOST_UBLAS_MATRIX_PROD_   
