AC_DEFUN([AX_GROMACS_LIBS], [
  AC_CHECK_LIB(gmx,GromacsVersion,[
    LIBS_GMX="-lgmx"
  ],[
    AC_CHECK_LIB(gmx_mpi,GromacsVersion,[
      LIBS_GMX="-lgmx_mpi -lmpi -llam -lpthread"
    ],[
      AC_MSG_ERROR([

Gromacs library not found!
please make sure LDFLAGS is pointing to <gromacs-path>/lib

])],[-lmpi -llam -lm] )],[-lm])
  AC_SUBST([LIBS_GMX])
])

AC_DEFUN([AX_GROMACS_HEADERS], [
AC_CHECK_HEADERS([copyrite.h],,[AC_MSG_ERROR([

Gromacs headers not found,
please make sure CPPFLAGS is pointing to <gomacs-path>/include])])

])

AC_DEFUN([AX_GROMACS], [
 AX_GROMACS_HEADERS
 AX_GROMACS_LIBS
])
