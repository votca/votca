AC_INIT(votca-tools,1.0_rc4,bugs@votca.org)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign])


AC_PROG_CXX
LT_INIT

AC_LANG_CPLUSPLUS

AC_HEADER_STDC
PKG_PROG_PKG_CONFIG

# check for boost
AX_BOOST_BASE([1.33.1])

# check for fftw3
PKG_CHECK_MODULES([FFTW],fftw3)

# check for libxml2
PKG_CHECK_MODULES([XML],libxml-2.0)

# check for gsl
PKG_CHECK_MODULES([GSL],gsl)

AC_CHECK_PROGS(TXT2TAGS,txt2tags,no)
AM_CONDITIONAL(HAVE_TXT2TAGS, test "$TXT2TAGS" != no)

AC_CONFIG_HEADERS(config.h)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libtools/Makefile])
AC_CONFIG_FILES([src/libtools/libvotca-tools.pc])

AC_OUTPUT
