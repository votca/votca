AC_DEFUN([AX_GROMACS_LIBS], [
#---------
  AC_CHECK_LIB(gmx,GromacsVersion,
    [LIBS_GMX="-lgmx"],[],[-lm])  

#----------
  if test -z "$LIBS_GMX"; then
    AC_CHECK_LIB(gmx_mpi,GromacsVersion,
      [LIBS_GMX="-lgmx_mpi -lmpi -lm -lpthread"],,[-lmpi -lpthread -lm])
  fi
#---------- 
  if test -z "$LIBS_GMX"; then
    AC_CHECK_LIB(gmx_mpi,GromacsVersion,
      [LIBS_GMX="-lgmx_mpi -lmpi -llam -lm -lpthread"],,[-lmpi -llam -lpthread -lm])
  fi

#------------
  if test -z "$LIBS_GMX"; then
    AC_MSG_ERROR([

Could not linka against GROMACS,
please check you LDFLAGS

    ])
   fi
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
