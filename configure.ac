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
            [AC_HELP_STRING([--with-votca_boost],
                            [Build libvotca_boost, a build-in replacement for boost])],,
                             [with_votca_boost="no"])
# check for boost
if test "$with_votca_boost" = no; then
  AX_BOOST_BASE([1.33.1],[have_boost=yes],[have_boost=no])
fi

if test "$have_boost" = "no" -o "$want_boost" = "no"; then
  AC_MSG_NOTICE([
		
BOOST was not found/disabled: Building libvotca_boost
])
  with_votca_boost="yes"
fi

# check for fftw3
PKG_CHECK_MODULES([FFTW],fftw3)

# check for libxml2
PKG_CHECK_MODULES([XML],libxml-2.0)

# check for gsl
PKG_CHECK_MODULES([GSL],gsl)

AC_CHECK_PROGS(TXT2TAGS,txt2tags,no)
AM_CONDITIONAL(HAVE_TXT2TAGS,[test "$TXT2TAGS" != no])

if test "$with_votca_boost" = "yes"; then
  BOOST_CPPFLAGS='-I$(top_srcdir)/src/libboost'
  AC_SUBST(BOOST_CPPFLAGS)
  BOOST_LDFLAGS='$(top_srcdir)/src/libboost/libs/libvotca_boost.la'
  AC_SUBST(BOOST_LDFLAGS)
  AC_CONFIG_SUBDIRS(src/libboost)
fi
AM_CONDITIONAL(VOTCA_BOOST, test "$with_votca_boost" = "yes")

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libtools/Makefile])
AC_CONFIG_FILES([src/libtools/libvotca_tools.pc])

AC_OUTPUT
