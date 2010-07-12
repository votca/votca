AC_INIT(votca-tools,1.0_rc5,bugs@votca.org)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.8 foreign])

AC_PROG_CXX
LT_INIT

AC_LANG([C++])

AC_HEADER_STDC
PKG_PROG_PKG_CONFIG

AC_CONFIG_HEADERS(config.h)

AC_ARG_ENABLE(votca-boost,
            [AS_HELP_STRING([--enable-votca-boost],
                            [Build libvotca_boost, a build-in replacement for boost])],,
                             [enable_votca_boost="no"])

AC_ARG_ENABLE(votca-expat,
            [AS_HELP_STRING([--enable-votca-expat],
                            [Build libvotca_expat, a build-in replacement for expat])],,
                             [enable_votca_expat="no"])

AC_ARG_ENABLE(la-files,
	[AS_HELP_STRING([--disable-la-files],[Do NOT install GNU libtool's la files])],,
	[enable_la_files=yes])

AM_CONDITIONAL(NO_LA_FILES,[test "$enable_la_files" = no])

# check for boost
if test "$enable_votca_boost" = no; then
  AX_BOOST_BASE([1.33.1],[have_boost=yes],[have_boost=no])
  if test "$have_boost" = "yes" -a "$want_boost" = "yes"; then
    AX_BOOST_PROGRAM_OPTIONS([1.33.1],[:],[have_boost=no])
  fi
  if test "$have_boost" = "no" -o "$want_boost" = "no"; then
    AC_MSG_ERROR([
		
BOOST was not found/disabled: Install BOOST or enable build-in replacement libvotca_boost with --enable-votca-boost
])
  fi
fi

# check for fftw3
PKG_CHECK_MODULES([FFTW],fftw3)

# check for gsl
PKG_CHECK_MODULES([GSL],gsl)

AC_CHECK_PROGS(TXT2TAGS,txt2tags,no)
AM_CONDITIONAL(HAVE_TXT2TAGS,[test "$TXT2TAGS" != no])

AC_CHECK_PROGS(DOXYGEN,doxygen,no)
AM_CONDITIONAL(HAVE_DOXYGEN,[test "$TXT2TAGS" != no])

if test "$enable_votca_expat" = "yes"; then
  EXPAT_CFLAGS='-I$(top_srcdir)/src/libexpat/include'
  EXPAT_LIBS='$(top_srcdir)/src/libexpat/libs/libvotca_expat.la'
  VOTCA_EXPAT="libexpat"
  AC_CONFIG_SUBDIRS(src/libexpat)
  PKGEXPAT="libvotca_expat"
  PKGCFLAGSEXPAT=""
  PKGLIBSEXPAT=""
  AC_SUBST(VOTCA_EXPAT)
  AC_SUBST(EXPAT_CFLAGS)
  AC_SUBST(EXPAT_LIBS)
  AC_SUBST(PKGCFLAGSEXPAT)
  AC_SUBST(PKGLIBSEXPAT)
  AC_SUBST(PKGEXPAT)
else
# check for expat
  AX_EXPAT
fi

if test "$enable_votca_boost" = "yes"; then
  BOOST_CFLAGS='-I$(top_srcdir)/src/libboost'
  BOOST_LIBS='$(top_srcdir)/src/libboost/libs/libvotca_boost.la'
  VOTCA_BOOST="libboost"
  AC_CONFIG_SUBDIRS(src/libboost)
  PKGBOOST="libvotca_boost"
  PKGCFLAGSBOOST=""
  PKGLIBSBOOST=""
else
  BOOST_CFLAGS="$BOOST_CPPFLAGS"
  BOOST_LIBS="$BOOST_LDFLAGS $BOOST_PROGRAM_OPTIONS_LIB"
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
AC_CONFIG_FILES([scripts/Makefile])
AC_CONFIG_FILES([share/Makefile])
AC_CONFIG_FILES([share/doc/Makefile])

AC_OUTPUT
