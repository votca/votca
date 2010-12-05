AC_DEFUN([AX_TXT2TAGS],[
  AC_ARG_WITH([txt2tags],[
    AS_HELP_STRING([--without-txt2tags@<:@=ARG@:>@],
      [disable usage of txt2tags needed for documentation and man pages.])
    ],,[with_txt2tags=yes])

  if test "$with_txt2tags" = "no"; then
    TXT2TAGS="no"
  else
    AC_CHECK_PROGS(TXT2TAGS,txt2tags,no)
    if test $TXT2TAGS = "no"; then
      AC_MSG_WARN([txt2tags not found, help configure to find it by setting TXT2TAGS])
    else
      WANT_TXT2TAGS_VERSION=ifelse([$1], ,2.5,$1)
      AC_MSG_CHECKING([for txt2tags version >= $WANT_TXT2TAGS_VERSION])
      HAVE_TXT2TAGS_VERSION=`$TXT2TAGS --version 2>&1 | sed 's/^.* \([[0-9.]]*\) .*$/\1/' 2> /dev/null`
      if expr "$WANT_TXT2TAGS_VERSION" \<= "$HAVE_TXT2TAGS_VERSION" > /dev/null 2>&1; then
        AC_MSG_RESULT([yes ($HAVE_TXT2TAGS_VERSION)])
        AC_MSG_CHECKING([if txt2tags works])
        if echo | $TXT2TAGS -t man -i - -o - > /dev/null 2>&1; then
          AC_MSG_RESULT([yes])
        else
          AC_MSG_RESULT([no])
          TXT2TAGS="no"
        fi
      else
        AC_MSG_RESULT([no ($HAVE_TXT2TAGS_VERSION)])
        TXT2TAGS="no"
      fi
    fi
  fi
  AM_CONDITIONAL(HAVE_TXT2TAGS,[test "$TXT2TAGS" != no])
])
