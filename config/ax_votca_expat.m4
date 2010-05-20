AC_DEFUN([AX_VOTCA_EXPAT], [
  PKG_CHECK_MODULES([VOTCA_EXPAT],libvotca_expat,[
    dnl we have votca_expat pkg-config file or variables
  ],[
    dnl we do not have votca_expat pkg-config file, try to VOTCALDLIB
    AC_ARG_VAR([VOTCALDLIB],[path to gromacs lib dir, usually set by "source VOTCARC"])
    AC_MSG_CHECKING([VOTCALDLIB])
    if test -z "$VOTCALDLIB"; then
      AC_MSG_RESULT([no])
    else
      AC_MSG_RESULT([yes])
      AC_MSG_NOTICE([creating VOTCA_EXPAT_LIBS and VOTCA_EXPAT_CFLAGS from VOTCALDLIB])
      if test -z "$VOTCA_EXPAT_LIBS"; then
        VOTCA_EXPAT_LIBS="-L$VOTCALDLIB -lvotca_expat"
        AC_MSG_NOTICE([setting VOTCA_EXPAT_LIBS   to "$VOTCA_EXPAT_LIBS"])
      else
        AC_MSG_NOTICE([VOTCA_EXPAT_LIBS was already set elsewhere to "$VOTCA_EXPAT_LIBS"])
      fi
      if test -z "$VOTCA_EXPAT_CFLAGS"; then
        VOTCA_EXPAT_CFLAGS="-I$VOTCALDLIB/../include/votca/expat"
        AC_MSG_NOTICE([setting VOTCA_EXPAT_CFLAGS to "$VOTCA_EXPAT_CFLAGS"])
      else
        AC_MSG_NOTICE([VOTCA_EXPAT_CFLAGS was already set elsewhere to "$VOTCA_EXPAT_CFLAGS"])
      fi
    fi
  ])
  save_CPPFLAGS="$CPPFLAGS"
  save_LIBS="$LIBS"

  CPPFLAGS="$VOTCA_EXPAT_CFLAGS $CPPFLAGS"
  LIBS="$VOTCA_EXPAT_LIBS $LIBS"
  AC_CHECK_HEADERS([expat.h],[
    AC_TRY_LINK_FUNC([XML_ParserCreate],[:],[ve_failed="yes"])
  ],[ve_failed="yes"])
  CPPFLAGS="$save_CPPFLAGS"
  LIBS="$save_LIBS"
  if test "$ve_failed" != "yes"; then
    PKG_CHECK_EXISTS(libvotca_expat,[
      PKGEXPAT="libvotca_expat"
      PKGCFLAGSEXPAT=""
      PKGLIBSEXPAT=""
    ],[
      PKGEXPAT=""
      PKGCFLAGSEXPAT="$VOTCA_EXPAT_CFLAGS"
      PKGLIBSEXPAT="$VOTCA_EXPAT_LIBS"
    ])
    EXPAT_CFLAGS="$VOTCA_EXPAT_CFLAGS"
    EXPAT_LIBS="$VOTCA_EXPAT_LIBS"
    AC_SUBST(EXPAT_CFLAGS)
    AC_SUBST(EXPAT_LIBS)
    AC_SUBST(PKGEXPAT)
    AC_SUBST(PKGLIBSEXPAT)
    AC_SUBST(PKGCFLAGSEXPAT)
  else
    dnl header or libs was missing, try system expat
    AX_EXPAT
  fi
])
