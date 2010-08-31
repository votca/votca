AC_DEFUN([AX_VOTCA_TOOLS], [
  PKG_CHECK_MODULES([VOTCA_TOOLS],libvotca_tools,[:],[
    AC_ARG_VAR([VOTCALDLIB],[path to gromacs lib dir, usually set by "source VOTCARC"])
    AC_MSG_CHECKING([VOTCALDLIB])
    if test -z "$VOTCALDLIB"; then
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([

Could not find VOTCALDLIB environment variable, please source <votca-path>/bin/VOTCARC.{csh,bash} 
or specify VOTCA_TOOLS_LIBS and VOTCA_TOOLS_CLFAGS
      ])
    else
      AC_MSG_RESULT([yes])
    fi
    AC_MSG_NOTICE([creating VOTCA_TOOLS_LIBS and VOTCA_TOOLS_CFLAGS from VOTCALDLIB])
    if test -z "$VOTCA_TOOLS_LIBS"; then
      VOTCA_TOOLS_LIBS="-L$VOTCALDLIB -lvotca_tools"
      AC_MSG_NOTICE([setting VOTCA_TOOLS_LIBS   to "$VOTCA_TOOLS_LIBS"])
    else
      AC_MSG_NOTICE([VOTCA_TOOLS_LIBS was already set elsewhere to "$VOTCA_TOOLS_LIBS"])
    fi
    if test -z "$VOTCA_TOOLS_CFLAGS"; then
      VOTCA_TOOLS_CFLAGS="-I$VOTCALDLIB/../include"
      AC_MSG_NOTICE([setting VOTCA_TOOLS_CFLAGS to "$VOTCA_TOOLS_CFLAGS"])
    else
      AC_MSG_NOTICE([VOTCA_TOOLS_CFLAGS was already set elsewhere to "$VOTCA_TOOLS_CFLAGS"])
    fi
  ])
  save_CPPFLAGS="$CPPFLAGS"
  save_LIBS="$LIBS"
  save_CXX="$CXX"

  CPPFLAGS="$VOTCA_TOOLS_CFLAGS $CPPFLAGS"
  LIBS="$VOTCA_TOOLS_LIBS $LIBS"
  AC_MSG_CHECKING([for votca::tools::ToolsVersionStr in $VOTCA_TOOLS_LIBS])
  CXX="${SHELL-/bin/sh} ./libtool --mode=link $CXX"
  AC_LINK_IFELSE([
    AC_LANG_PROGRAM([#include <votca/tools/version.h>],[votca::tools::ToolsVersionStr()])
  ],[
    AC_MSG_RESULT([yes])
  ],[
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([

Could not link against VOTCA tools,
please check your LDFLAGS and/or specify libraries required to link 
in VOTCA_TOOLS_LIBS (e.g. export VOTCA_TOOLS_LIBS="-L<votca-path>/lib -lvotca_tools").
    ])
  ])
  CXX="$save_CXX"

  AC_CHECK_HEADERS([votca/tools/version.h],,[
    AC_MSG_ERROR([

Votca tools headers not found,
please make sure VOTCA_TOOLS_CFLAGS is pointing to <votca-path>/include
(e.g. export VOTCA_TOOLS_CFLAGS="-I<votca-path>/include").
    ])
  ])

  AC_CHECK_HEADERS([votca/tools/application.h],,[
    AC_MSG_ERROR([

Votca tools headers were found, but boost headers not!
please make sure that VOTCA_TOOLS_CFLAGS is pointing to the votca tools headers AND to the boost headers
(e.g. export VOTCA_TOOLS_CFLAGS="-I<votca-path>/include -I<path/to/boost>/include").

If you are using votca-boost (build-in replacement for boost) please export
VOTCA_TOOLS_CFLAGS="-I<votca-path>/include -I<votca-path>/include/votca".
    ])
  ])
  CPPFLAGS="$save_CPPFLAGS"
  LIBS="$save_LIBS"
])
