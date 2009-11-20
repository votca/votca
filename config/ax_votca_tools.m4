AC_DEFUN([AX_VOTCA_TOOLS_LIBS], [
#---------
  AC_TRY_LINK([#include <tools/property.h>], [ Property p; return 0; ], [
      LIBS_TOOLS="-ltools" ],[])

#AX_TRY_LIB([tools],[load_property_from_xml],
#    LIBS_TOOLS="-ltools",,[$LIBS_GSL])  


#------------
  if test -z "$LIBS_TOOLS"; then
    AC_MSG_ERROR([

Could not link against votca tools library,
please make sure you have it installed and check your LDFLAGS

    ])
   fi
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
