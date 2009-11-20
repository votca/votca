AC_DEFUN([AX_TRY_LIB],[
  	try_lib="$1"
  	try_func="$2"
        pushdef([try_if_yes], [$3])
	pushdef([try_if_no], [$4])
	try_deps="$5"
	

  try_libs_save="$LIBS"
  LIBS="-l$try_lib $try_deps $LIBS"

  AC_MSG_CHECKING([for $try_func with -l$try_lib $try_deps"])
  AC_TRY_LINK_FUNC($try_func,try_found="yes",try_found="no")
  AC_MSG_RESULT($try_found)

  LIBS="$try_libs_save"

  if test "x$try_found" = "xno" ; then
    ifelse(try_if_no, , echo >/dev/null, try_if_no)
  else
    ifelse(try_if_yes, , echo >/dev/null, try_if_yes)
  fi
])
