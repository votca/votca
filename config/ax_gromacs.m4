AC_DEFUN([AX_GROMACS], [
  PKG_CHECK_MODULES([GMX],libgmx)
  save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$GMX_CFLAGS $CPPFLAGS"
  AC_CHECK_HEADERS([tpxio.h],,[AC_MSG_ERROR([

Gromacs headers not found,
please make sure GMX_CFLAGS is pointing to <gomacs-path>/include/gromacs 
Do not forget the /gromacs due to bug in gromacs headers!])])
 
 CPPFLAGS="$save_CPPFLAGS"
 AC_CHECK_LIB([gmx],GromacsVersion,[HAVE_GMX_LIBS="yes"],[HAVE_GMX_LIBS="no"],[$GMX_LIBS])
 if test "HAVE_GMX_LIBS" = "no"; then
   AC_CHECK_LIB([gmx_mpi],GromacsVersion,[HAVE_GMX_LIBS="yes"],[HAVE_GMX_LIBS="no"],[$GMX_LIBS])
 fi
 if test "HAVE_GMX_LIBS" = "no"; then
   AC_MSG_ERROR([

Could not link against GROMACS,
please check your LDFLAGS and/or specify libraries required to link 
against GROMACS in GMX_LIBS (e.g. export GMX_LIBS="-Lpath/to/gramacs/lib -lgmx").

If you are using a mpi version of gromacs, make sure that CXX is something like mpic++.

   ])
 fi

 save_CPPFLAGS="$CPPFLAGS"
 CPPFLAGS="$GMX_CFLAGS $CPPFLAGS"
 AC_CHECK_HEADERS([types/oenv.h],[gmx_dev=yes],[gmx_dev=no])
 CPPFLAGS="$save_CPPFLAGS"
 if test "$gmx_dev" = "yes"; then
   AC_CHECK_LIB([gmx],output_env_init,[gmx_dev="yes"],[AC_MSG_ERROR([

We found oenv.h header, which exist only in gromacs 4 devel, but we did find
output_env_init libgmx

])],[$GMX_LIBS])
   AC_MSG_NOTICE([

We are using a development version of gromacs, we hope you know what you are doing....
unexpected results ahead

])
 AC_DEFINE(GMX4DEV,,[Use gromacs 4 devel version])
 fi
])
