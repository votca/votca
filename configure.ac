AC_INIT(votca-tools,1.0-dev,bugs@votca.org)

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

AC_ARG_WITH(fftw,
            [AS_HELP_STRING([--without-fftw],
                            [Compile without fftw, will disable CrossCorrelate class @<:@ARG=no@:>@])],,
                             [with_fftw="yes"])

AC_ARG_WITH(gsl,
            [AS_HELP_STRING([--without-gsl],
                            [Compile without gsl, will disable linalg functions and CubicSpline class @<:@ARG=no@:>@])],,
                             [with_gsl="yes"])

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

if test "$with_fftw" = "yes"; then
  # check for fftw3
  PKG_CHECK_MODULES([FFTW],fftw3)
  AC_SUBST(PKGFFTW,"fftw3")
else
  AC_DEFINE([NOFFTW],,[Compile without fftw and disable CrossCorrelate class])
  AC_SUBST(FFTW_CFLAGS,"")
  AC_SUBST(FFTW_LIBS,"")
  AC_SUBST(PKGFFTW,"")
fi

AC_CHECK_PROGS(TXT2TAGS,txt2tags,no)
AM_CONDITIONAL(HAVE_TXT2TAGS,[test "$TXT2TAGS" != no])

AC_CHECK_PROGS(DOXYGEN,doxygen,no)
AM_CONDITIONAL(HAVE_DOXYGEN,[test "$DOXYGEN" != no])

if test "$enable_votca_expat" = "yes"; then
  AC_CONFIG_SUBDIRS(src/libexpat)
  AC_SUBST(VOTCA_EXPAT,"libexpat")
  AC_SUBST(EXPAT_CFLAGS,'-I$(top_srcdir)/src/libexpat/include')
  AC_SUBST(EXPAT_LIBS,'$(top_builddir)/src/libexpat/libs/libvotca_expat.la')
  AC_SUBST(PKGEXPAT,"libvotca_expat")
  AC_SUBST(PKGCFLAGSEXPAT,"")
  AC_SUBST(PKGLIBSEXPAT,"")
else
  AC_SUBST(VOTCA_EXPAT,"")
  # check for expat
  AX_EXPAT
fi

if test "$with_gsl" = "yes"; then
  # check for gsl
  PKG_CHECK_MODULES([GSL],gsl)
  AC_SUBST(PKGGSL,"gsl")
else
  AC_DEFINE([NOGSL],,[Compile without gsl and disable inalg functions and CubicSpline class])
  AC_SUBST(GSL_CFLAGS,"")
  AC_SUBST(GSL_LIBS,"")
  AC_SUBST(PKGGSL,"")
fi

if test "$enable_votca_boost" = "yes"; then
  AC_ARG_VAR([VOTCA_BOOST_CXXFLAGS],[C++ compiler flags for libvotca_boost, overwrites its CXXFLAGS])
  AC_CONFIG_SUBDIRS(src/libboost)
  AC_SUBST(VOTCA_BOOST,"libboost")
  AC_SUBST(BOOST_CFLAGS,'-I$(top_srcdir)/src/libboost')
  AC_SUBST(BOOST_LIBS,'$(top_builddir)/src/libboost/libs/libvotca_boost.la')
  AC_SUBST(PKGBOOST,"libvotca_boost")
  AC_SUBST(PKGCFLAGSBOOST,"")
  AC_SUBST(PKGLIBSBOOST,"")
else
  AC_SUBST(VOTCA_BOOST,"")
  AC_SUBST(BOOST_CFLAGS,"$BOOST_CPPFLAGS")
  AC_SUBST(BOOST_LIBS,"$BOOST_LDFLAGS $BOOST_PROGRAM_OPTIONS_LIB")
  AC_SUBST(PKGBOOST,"")
  AC_SUBST(PKGCFLAGSBOOST,"$BOOST_CFLAGS")
  AC_SUBST(PKGLIBSBOOST,"$BOOST_LIBS")
fi

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libtools/Makefile])
AC_CONFIG_FILES([src/libtools/libvotca_tools.pc])
AC_CONFIG_FILES([scripts/Makefile])
AC_CONFIG_FILES([share/Makefile])
AC_CONFIG_FILES([share/doc/Makefile])

AC_OUTPUT
