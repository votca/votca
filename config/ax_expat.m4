AC_DEFUN([AX_EXPAT], [
  dnl next version of expat have pkg-config support
  PKG_CHECK_MODULES([EXPAT],expat,[
    PKGEXPAT="expat"
    PKGCFLAGSEXPAT=""
    PKGLIBSEXPAT=""
  ],[
    save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$EXPAT_CFLAGS $CPPFLAGS"
    AC_CHECK_HEADER([expat.h],[:],[AC_MSG_ERROR([
Expat headers not found,
please make sure EXPAT_CFLAGS is pointing to <expat-path>/include])
    ]) 
    CPPFLAGS="$save_CPPFLAGS"
    if test -z "$EXPAT_LIBS"; then
      EXPAT_LIBS="-lexpat"
    fi
    save_LIBS="$LIBS"
    LIBS="$EXPAT_LIBS $LIBS"
    AC_MSG_CHECKING([for XML_ParserCreate in $EXPAT_LIBS])
    AC_TRY_LINK_FUNC([XML_ParserCreate],[AC_MSG_RESULT([yes])],[
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([

Could not link against Expat,
please check your LDFLAGS and/or specify libraries required to link 
against expat in EXPAT_LIBS (e.g. export EXPAT_LIBS="-L<expat-pat>/lib -lexpat").
      ])
    ])
    LIBS="$save_LIBS"
    PKGEXPAT="expat"
    PKGCFLAGSEXPAT="$EXPAT_CFLAGS"
    PKGLIBSEXPAT="$EXPAT_LIBS"
  ])
  AC_SUBST(PKGCFLAGSEXPAT)
  AC_SUBST(PKGLIBSEXPAT)
  AC_SUBST(PKGEXPAT)
])
