AC_INIT(moo,0.5,werauchimmer@mpip-mainz.mpg.de)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_HEADERS(config.h)

AC_PROG_CXX
LT_INIT

AC_LANG_CPLUSPLUS

AC_HEADER_STDC
PKG_PROG_PKG_CONFIG

PKG_CHECK_MODULES([VOTCA_TOOLS],libvotca_tools)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libmoo/Makefile])
AC_CONFIG_FILES([src/libmoo/libvotca_moo.pc])

AC_OUTPUT

