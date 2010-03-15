AC_DEFUN([AX_VOTCA_TOOLS], [
  PKG_CHECK_MODULES([VOTCA_TOOLS],libvotca_tools,[:],[
    AC_ARG_VAR([VOTCALDLIB],[path to gromacs lib dir, usually set by "source GMXRC"])
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

  CPPFLAGS="$VOTCA_TOOLS_CFLAGS $CPPFLAGS"
  LIBS="$VOTCA_TOOLS_LIBS $LIBS"
  AC_CHECK_HEADERS([votca/tools/property.h],,[
    AC_MSG_ERROR([

Votca tools headers not found,
please make sure VOTCA_TOOLS_CFLAGS is pointing to <votca-path>/include
    ])
  ])
  AC_MSG_CHECKING([for exit in $VOTCA_TOOLS_LIBS])
  AC_TRY_LINK_FUNC(exit,[AC_MSG_RESULT([yes])],[
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([

Could not link against VOTCA tools,
please check your LDFLAGS and/or specify libraries required to link 
in VOTCA_TOOLS_LIBS (e.g. export VOTCA_TOOLS_LIBS="-L<votca-path>/lib -lvotca_tools").
    ])
  ])
  CPPFLAGS="$save_CPPFLAGS"
  LIBS="$save_LIBS"
])
