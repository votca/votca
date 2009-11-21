AC_DEFUN([AX_VOTCA_TOOLS_LIBS], [
#---------
  tmp_libs="$LIBS"
  LIBS="-ltools $LIBS_XML $LIBS"
  AC_MSG_CHECKING([for votca tools library])
  AC_TRY_LINK([#include <tools/property.h>], [ Property p; p.get("foo"); ], [
      LIBS_TOOLS="-ltools" ],[AC_MSG_RESULT([no])])
  LIBS="$tmp_libs"

#------------
  if test -z "$LIBS_TOOLS"; then
    AC_MSG_ERROR([

Could not link against votca tools library,
please make sure you have it installed and check your LDFLAGS

    ])
   fi

   AC_MSG_RESULT([yes])
   AC_SUBST([LIBS_TOOLS])
])

AC_DEFUN([AX_VOTCA_TOOLS_HEADERS], [
  AC_CHECK_HEADERS([tools/property.h],,[
    AC_MSG_ERROR([

Header from votca tools library not found!
please make sure you have it installed and set your CPPFLAGS properly!

    ])
  ])

])

AC_DEFUN([AX_VOTCA_TOOLS], [
 AX_VOTCA_TOOLS_HEADERS
 AX_VOTCA_TOOLS_LIBS
])
