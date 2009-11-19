AC_DEFUN([AX_TRY_LIB],[
  	try_lib="$1"
  	try_func="$2"
        pushdef([try_if_yes], [$3])
	pushdef([try_if_no], [$4])
	try_deps="$5"
	

  try_libs_save="$LIBS"
  LIBS="-l$try_lib $try_deps $LIBS"

  AC_MSG_CHECKING([for $try_func with -l$try_lib $try_deps"])
  AC_TRY_LINK_FUNC($try_func,try_found="yes",try_found="no")
  AC_MSG_RESULT($try_found)

  LIBS="$try_libs_save"

  if test "x$try_found" = "xno" ; then
    ifelse(try_if_no, , echo >/dev/null, try_if_no)
  else
    ifelse(try_if_yes, , echo >/dev/null, try_if_yes)
  fi
])

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
AC_CHECK_HEADERS([copyrite.h],,[AC_MSG_ERROR([

Gromacs headers not found,
please make sure CPPFLAGS is pointing to <gomacs-path>/include])])

])

AC_DEFUN([AX_GROMACS], [
 AX_GROMACS_HEADERS
 AX_GROMACS_LIBS
])
