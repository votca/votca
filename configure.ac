AC_INIT(moo,0.5,werauchimmer@mpip-mainz.mpg.de)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_HEADERS(config.h)

AC_PROG_CXX
AC_PROG_RANLIB

AC_LANG_CPLUSPLUS

AC_HEADER_STDC
PKG_PROG_PKG_CONFIG

#check for boost
AX_BOOST_BASE([1.33.1])

if test "$have_boost" = "no" -o "$want_boost" = "no"; then
  AC_MSG_NOTICE([
		
BOOST was not found/disabled: Try to use libvotca_boost

])
  PKG_CHECK_MODULES([BOOST],libvotca_boost)
else
  BOOST_CFLAGS="$BOOST_CPPFLAGS"
  AC_SUBST(BOOST_CFLAGS)
  BOOST_LIBS="$BOOST_LDFLAGS"
  AC_SUBST(BOOST_LIBS)
fi

#look for votca tools
PKG_CHECK_MODULES([VOTCA_TOOLS],libvotca_tools)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libmoo/Makefile])
AC_CONFIG_FILES([src/libmoo/libvotca_moo.pc])
AC_CONFIG_FILES([src/tools/Makefile])

AC_OUTPUT

