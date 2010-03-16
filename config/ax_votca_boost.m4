AC_DEFUN([AX_VOTCA_BOOST], [
  PKG_CHECK_MODULES([VOTCA_BOOST],libvotca_boost,[
    dnl we have votca_boost pkg-config file or variables
  ],[
    dnl we do not have votca_boost pkg-config file, try to VOTCALDLIB
    AC_ARG_VAR([VOTCALDLIB],[path to gromacs lib dir, usually set by "source VOTCARC"])
    AC_MSG_CHECKING([VOTCALDLIB])
    if test -z "$VOTCALDLIB"; then
      AC_MSG_RESULT([no])
    else
      AC_MSG_RESULT([yes])
      AC_MSG_NOTICE([creating VOTCA_BOOST_LIBS and VOTCA_BOOST_CFLAGS from VOTCALDLIB])
      if test -z "$VOTCA_BOOST_LIBS"; then
        VOTCA_BOOST_LIBS="-L$VOTCALDLIB -lvotca_boost"
        AC_MSG_NOTICE([setting VOTCA_BOOST_LIBS   to "$VOTCA_BOOST_LIBS"])
      else
        AC_MSG_NOTICE([VOTCA_BOOST_LIBS was already set elsewhere to "$VOTCA_BOOST_LIBS"])
      fi
      if test -z "$VOTCA_BOOST_CFLAGS"; then
        VOTCA_BOOST_CFLAGS="-I$VOTCALDLIB/../include/votca"
        AC_MSG_NOTICE([setting VOTCA_BOOST_CFLAGS to "$VOTCA_BOOST_CFLAGS"])
      else
        AC_MSG_NOTICE([VOTCA_BOOST_CFLAGS was already set elsewhere to "$VOTCA_BOOST_CFLAGS"])
      fi
    fi
  ])
  save_CPPFLAGS="$CPPFLAGS"
  save_LIBS="$LIBS"

  CPPFLAGS="$VOTCA_BOOST_CFLAGS $CPPFLAGS"
  LIBS="$VOTCA_BOOST_LIBS $LIBS"
  AC_CHECK_HEADERS([boost/program_options.hpp],[
    AC_MSG_CHECKING([for boost::program_options::value in $VOTCA_BOOST_LIBS])
    dnl no libtool wrapper needed here, because boost has no deps itself
    AC_LINK_IFELSE([
      AC_LANG_PROGRAM([#include <boost/program_options.hpp>
      ],[boost::program_options::value<int>()])
    ],[
      AC_MSG_RESULT([yes])
    ],[
      AC_MSG_RESULT([no])
      vb_failed="yes"
    ])
  ],[vb_failed="yes"])
  CPPFLAGS="$save_CPPFLAGS"
  LIBS="$save_LIBS"
  if test "$vb_failed" != "yes"; then
    PKG_CHECK_EXISTS(libvotca_boost,[
      PKGBOOST="libvotca_boost"
      PKGCFLAGSBOOST=""
      PKGLIBSBOOST=""
    ],[
      PKGBOOST=""
      PKGCFLAGSBOOST="$VOTCA_BOOST_CFLAGS"
      PKGLIBSBOOST="$VOTCA_BOOST_LIBS"
    ])
    BOOST_CFLAGS="$VOTCA_BOOST_CFLAGS"
    BOOST_LIBS="$VOTCA_BOOST_LIBS"
  else
    dnl header or libs was missing, try system boost
    AX_BOOST_BASE([1.33.1],[:],[
      AC_MSG_ERROR([
  		
No system BOOST and no libvotca_boost found.
      ])
    ])
    if test "$want_boost" = "no"; then
      AC_MSG_ERROR([
  		
You can not disable boost without having libvotca_boost.
      ])
    fi
    AX_BOOST_PROGRAM_OPTIONS([1.33.1])
    BOOST_CFLAGS="$BOOST_CPPFLAGS"
    BOOST_LIBS="$BOOST_LDFLAGS $BOOST_PROGRAM_OPTIONS_LIB"
    PKGBOOST=""
    PKGLIBSBOOST="$BOOST_LIBS"
    PKGCFLAGSBOOST="$BOOST_CFLAGS"
  fi
  AC_SUBST(BOOST_CFLAGS)
  AC_SUBST(BOOST_LIBS)
  AC_SUBST(PKGBOOST)
  AC_SUBST(PKGLIBSBOOST)
  AC_SUBST(PKGCFLAGSBOOST)
])
