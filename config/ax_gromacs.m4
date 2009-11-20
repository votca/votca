AC_DEFUN([AX_GROMACS_LIBS], [
#---------
  AX_TRY_LIB([gmx],[GromacsVersion],
    LIBS_GMX="-lgmx",,[-lm])  

#----------
  if test -z "$LIBS_GMX"; then
    AX_TRY_LIB(gmx_mpi,GromacsVersion,
      [LIBS_GMX="-lgmx_mpi -lmpi -lm -lpthread"],[],[-lmpi -lpthread -lm])
  fi
#---------- 
  if test -z "$LIBS_GMX"; then
    AX_TRY_LIB(gmx_mpi,GromacsVersion,
      [LIBS_GMX="-lgmx_mpi -lmpi -llam -lm -lpthread"],,[-lmpi -llam -lpthread -lm])
  fi

#------------
  if test -z "$LIBS_GMX"; then
    AC_MSG_ERROR([

Could not link against GROMACS,
please check you LDFLAGS

    ])
   fi
   AC_SUBST([LIBS_GMX])
])

AC_DEFUN([AX_GROMACS_HEADERS], [
AC_CHECK_HEADERS([tpxio.h],,[AC_MSG_ERROR([

Gromacs headers not found,
please make sure CPPFLAGS is pointing to <gomacs-path>/include/gromacs 
Do not forget the /gromacs due to bug in gromacs headers!])])

])

AC_DEFUN([AX_GROMACS], [
 AX_GROMACS_HEADERS
 AX_GROMACS_LIBS
])
