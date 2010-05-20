AC_INIT(votca-expat,2.0.1,bugs@votca.org)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.8 foreign])

AC_PROG_CC
LT_INIT

AC_LANG([C])

AC_HEADER_STDC

AC_CONFIG_HEADER(expat_config.h)

AC_C_CONST
AC_TYPE_SIZE_T
AC_CHECK_FUNCS(memmove bcopy)

dnl Some basic configuration:
AC_DEFINE([XML_NS], 1,
          [Define to make XML Namespaces functionality available.])
AC_DEFINE([XML_DTD], 1,
          [Define to make parameter entity parsing functionality available.])
AC_DEFINE([XML_CONTEXT_BYTES], 1024,
          [Define to specify how much context to retain around the current parse point.])

AC_ARG_ENABLE(la-files,
	[AS_HELP_STRING([--disable-la-files],[Do NOT install GNU libtool's la files])],,
	[enable_la_files=yes])

AM_CONDITIONAL(NO_LA_FILES,[test "$enable_la_files" = no])

AC_CONFIG_FILES([libvotca_expat.pc])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([libs/Makefile])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
