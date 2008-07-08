# copied from GROMACS (www.gromacs.org)
# Much simpler check than fftw2
# Check for header <fftw3.h> AC_CHECK_HEADERS doesnt work, since we must
# use mpicc to get includes - cpp isnt always the same compiler.

AC_DEFUN([AX_FFTW3], [
  AC_MSG_CHECKING([for fftw3.h])
  AC_TRY_COMPILE([#include<fftw3.h>],,
  [
    # ok, look for library file too
    AC_MSG_RESULT(yes)
    if test "$enable_float" = "yes"; then
      AC_CHECK_LIB([fftw3f],main,,AC_MSG_ERROR([Cannot find fftw3f library]))
    else
      AC_CHECK_LIB([fftw3],main,,AC_MSG_ERROR([Cannot find fftw3 library]))
    fi
      MDLIB_LIBOBJS="$MDLIB_LIBOBJS gmx_fft_fftw3.lo"
  ],
  [
    # not ok, echo a warning
    AC_MSG_ERROR(
    [Cannot find FFT library (fftw3).

      Use CPPFLAGS and LDFLAGS if the library is installed in a
      non-standard location. ])
  ])
])

