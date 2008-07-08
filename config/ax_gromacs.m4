AC_DEFUN([AX_GROMACS_LIBS], [
AC_CHECK_LIB(gmx,GromacsVersion,,[AC_MSG_ERROR([

Gromacs library not found!
please make sure LDFLAGS is pointing to <gromacs-path>/lib

])])
])

AC_DEFUN([AX_GROMACS_HEADERS], [
AC_CHECK_HEADERS([gromacs/copyrite.h],,[AC_MSG_ERROR([

Gromacs headers not found,
please make sure CPPFLAGS is pointing to <gomacs-path>/include])])

])

AC_DEFUN([AX_GROMACS], [
 AX_GROMACS_HEADERS
 AX_GROMACS_LIBS
])
