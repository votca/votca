AC_INIT(moo,0.5,werauchimmer@mpip-mainz.mpg.de)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_HEADERS(config.h)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/libmoo/Makefile])

AC_PROG_CXX
AC_PROG_RANLIB

AC_LANG_CPLUSPLUS

AC_HEADER_STDC

AC_PROG_MAKE_SET

#check for boost
AX_BOOST([1.33.1])

ACI_PACKAGE_SAVED(XML,[XML library (version XML2)],
                  xml2-config,,,
                  AC_MSG_WARN(Cannot find the XML2 libs. Is xml2-config in your path ?)
                 )

#AX_BOOST_PROGRAM_OPTIONS
AX_VOTCA_TOOLS

AC_OUTPUT

