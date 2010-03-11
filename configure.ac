AC_INIT(votca-tools,1.0_rc4,bugs@votca.org)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.8 foreign])


AC_PROG_CXX
LT_INIT

AC_LANG_CPLUSPLUS

AC_HEADER_STDC
PKG_PROG_PKG_CONFIG

AC_CONFIG_HEADERS(config.h)

AC_ARG_WITH(votca_boost,
            [AC_HELP_STRING([--with-votca-boost],
                            [Build libvotca_boost, a build-in replacement for boost])],,
                             [with_votca_boost="no"])
# check for boost
if test "$with_votca_boost" = no; then
  AX_BOOST_BASE([1.33.1],[:],[have_boost=no])
  if test "$have_boost" = "no" -o "$want_boost" = "no"; then
    AC_MSG_NOTICE([
		
BOOST was not found/disabled: Building libvotca_boost
])
    with_votca_boost="yes"
  fi
fi

# check for fftw3
PKG_CHECK_MODULES([FFTW],fftw3)

# checko for expat
AX_EXPAT

# check for gsl
PKG_CHECK_MODULES([GSL],gsl)

AC_CHECK_PROGS(TXT2TAGS,txt2tags,no)
AM_CONDITIONAL(HAVE_TXT2TAGS,[test "$TXT2TAGS" != no])

if test "$with_votca_boost" = "yes"; then
  BOOST_CFLAGS='-I$(top_srcdir)/src/libboost'
  BOOST_LIBS='$(top_srcdir)/src/libboost/libs/libvotca_boost.la'
  VOTCA_BOOST="libboost"
  AC_CONFIG_SUBDIRS(src/libboost)
  PKGBOOST="libvotca_boost"
  PKGCFLAGSBOOST=""
  PKGLIBSBOOST=""
else
  BOOST_CFLAGS="$BOOST_CPPFLAGS"
  BOOST_LIBS=""
  VOTCA_BOOST=""
  PKGBOOST=""
  PKGCFLAGSBOOST="$BOOST_CFLAGS"
  PKGLIBSBOOST="$BOOST_LIBS"
fi
AC_SUBST(VOTCA_BOOST)
AC_SUBST(BOOST_CFLAGS)
AC_SUBST(BOOST_LIBS)
AC_SUBST(PKGCFLAGSBOOST)
AC_SUBST(PKGLIBSBOOST)
AC_SUBST(PKGBOOST)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libtools/Makefile])
AC_CONFIG_FILES([src/libtools/libvotca_tools.pc])

AC_OUTPUT
