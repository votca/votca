AC_DEFUN([AX_VOTCA_BOOST], [
  PKG_CHECK_MODULES([VOTCA_BOOST],libvotca_boost,[
    dnl we have votca_boost		   
    BOOST_CFLAGS="$VOTCA_BOOST_CFLAGS"
    BOOST_LIBS="$VOTCA_BOOST_LIBS"
    PKGBOOST="libvotca_boost"
    PKGLIBSBOOST=""
    PKGCFLAGSBOOST=""
  ],[
    dnl we do not have votca_boost, try system boost
    AX_BOOST_BASE([1.33.1],[:],[
      AC_MSG_ERROR([
  		
No system BOOST and no libvotca_boost found.
      ])
    ])
    if test "$want_boost" = "no"; then
      AC_MSG_ERROR([
  		
You can not disable boost without having libvotca_boost.
      ])
    fi
    AX_BOOST_PROGRAM_OPTIONS([1.33.1])
    BOOST_CFLAGS="$BOOST_CPPFLAGS"
    BOOST_LIBS="$BOOST_LDFLAGS $BOOST_PROGRAM_OPTIONS_LIB"
    PKGBOOST=""
    PKGLIBSBOOST="$BOOST_PROGRAM_OPTIONS_LIB"
    PKGCFLAGSBOOST="$BOOST_PROGRAM_OPTIONS_LIB"
  ])
  AC_SUBST(BOOST_CFLAGS)
  AC_SUBST(BOOST_LIBS)
  AC_SUBST(PKGBOOST)
  AC_SUBST(PKGLIBSBOOST)
  AC_SUBST(PKGCFLAGSBOOST)
])
