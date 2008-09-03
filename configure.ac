AC_INIT(csg,0.5,ruehle@mpip-mainz.mpg.de)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)

AM_INIT_AUTOMAKE([1.9 foreign])

AC_CONFIG_HEADERS(src/config.h)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([src/csg/Makefile])
AC_CONFIG_FILES([include/Makefile])
AC_CONFIG_FILES([src/tools/Makefile])
AC_CONFIG_FILES([src/libcsg/Makefile])

AC_PROG_CXX
AC_PROG_RANLIB
AC_HEADER_STDC
AC_PROG_MAKE_SET

#check for boost
AX_BOOST([1.33.1])
#AX_BOOST_PROGRAM_OPTIONS

# check for fftw3
AX_FFTW3

# check for libxml2
ACI_PACKAGE_SAVED(XML,[XML library (version XML2)],
                  xml2-config,,,
                  AC_MSG_WARN(Cannot find the XML2 libs. Is xml2-config in your path ?)
                 )

# check for gromacs
# AX_GROMACS_HEADERS

AC_OUTPUT
