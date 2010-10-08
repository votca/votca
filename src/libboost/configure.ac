AC_INIT(votca-boost,1.42.0-0,bugs@votca.org)

AC_CONFIG_AUX_DIR(../../config)
AC_CONFIG_MACRO_DIR(../../config)

AM_INIT_AUTOMAKE([1.8 foreign])

AC_PROG_CXX
LT_INIT

AC_LANG([C++])

AC_HEADER_STDC

AC_CONFIG_HEADERS(config.h)

AC_ARG_ENABLE(la-files,
	[AS_HELP_STRING([--disable-la-files],[Do NOT install GNU libtool's la files])],,
	[enable_la_files=yes])

AM_CONDITIONAL(NO_LA_FILES,[test "$enable_la_files" = no])

AC_CONFIG_FILES([libvotca_boost.pc])
AC_CONFIG_FILES([boost/Makefile])
AC_CONFIG_FILES([libs/Makefile])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
