AC_DEFUN([AX_GROMACS_LIBS], [
#---------
  HAVE_GMX_LIBS="no"

  if test ! -z "$LIBS_GMX"; then
    AX_TRY_LIB(,[GromacsVersion],
    HAVE_GMX_LIBS="yes",,[$LIBS_GMX])
  fi
#---------
  if test -z "$LIBS_GMX"; then
    AX_TRY_LIB([gmx],[GromacsVersion],[
      LIBS_GMX="-lgmx"
      HAVE_GMX_LIBS="yes"],,[-lm])  
  fi

#----------
  if test -z "$LIBS_GMX"; then
    AX_TRY_LIB(gmx_mpi,GromacsVersion,[
      LIBS_GMX="-lgmx_mpi -lm";
      HAVE_GMX_LIBS="yes"],,[-lm])
  fi
#------------
  if test "x$HAVE_GMX_LIBS" = "xno"; then
    AC_MSG_ERROR([

Could not link against GROMACS,
please check your LDFLAGS and/or specify libraries required to link 
against GROMACS in LIBS_GMX (e.g. export LIBS_GMX="-lgmx -lpthread -lxml2 -lm").

If you are using a mpi version of gromacs, make sure that CXX is something like mpic++.

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
