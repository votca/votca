#ifndef __VOTCA_TOOLS_VOTCA_CONFIG_H
#define __VOTCA_TOOLS_VOTCA_CONFIG_H

/* Compile without fftw and disable CrossCorrelate class */
/* #undef NOFFTW */

/* OpenMP */

#if defined(_OPENMP)
   #include <omp.h>
#endif



/* Linear algebra packages */
/* #undef EIGEN */
#define MKL
/* #undef GSL */




#if defined(GSL)
    #include "gsl_boost_ublas_matrix_prod.h"
    #ifndef NDEBUG
    #define NDEBUG
    #endif
#elif defined(MKL)
    #include "mkl_boost_ublas_matrix_prod.h"
    #ifndef NDEBUG
    #define NDEBUG
    #endif
#endif

/* Version number of package */
#define TOOLS_VERSION "1.4-dev"

#endif // __VOTCA_TOOLS_VOTCA_CONFIG_H
