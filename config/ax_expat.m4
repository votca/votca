AC_DEFUN([AX_EXPAT], [
  dnl next version of expat have pkg-config support
  PKG_CHECK_MODULES([EXPAT],expat,[:],[
    save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$EXPAT_CFLAGS $CPPFLAGS"
    AC_CHECK_HEADERS([expat.h],[:],[AC_MSG_ERROR([
Expat headers not found,
please make sure EXPAT_CFLAGS is pointing to <expat-path>/include])
    ]) 
    CPPFLAGS="$save_CPPFLAGS"
    if test -z "$EXPAT_LIBS"; then
      EXPAT_LIBS="-lexpat"
    fi
    AC_CHECK_LIB([expat],XML_ParserCreate,[:],[AC_MSG_ERROR([

Could not link against Expat,
please check your LDFLAGS and/or specify libraries required to link 
against expat in EXPAT_LIBS (e.g. export EXPAT_LIBS="-Lpath/to/gramacs/lib -lexpat").

If you are using a mpi version of gromacs, make sure that CXX is something like mpic++.

   ])],[$EPXAT_LIBS])
  ])
])
