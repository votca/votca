#! /bin/bash

die() {
  echo -e "$*" >&2
  exit 1
}

for i in *; do
  [ "$i" = "libs" ] && continue
  [ "$i" = "boost" ] && continue
  [ -d "$i" ] && die "Unknown dir $i"
done

echo Doing libs
cd libs
trunc_dir="$PWD"
files=$(find . -type f -not -name "*.cpp" -a -not -name "Makefile*")
[ -n "$files" ] && die "Unknown files:\n$files"

for i in $(find . -type d); do
  cd "$trunc_dir"
  echo Entering $i
  cd $i
  [ "$(echo *)" = '*' ] && continue
  rm -f Makefile.am
  dirs=""
  for j in *; do
    [ -d "$j" ] && dirs="$dirs $j" 
  done
  [ -n "$dirs" ] && echo -e "SUBDIRS = $dirs\n" >> Makefile.am
  [ "$(echo *.cpp)" = '*.cpp' ] && continue
  echo -e "EXTRA_DIST = \\" >> Makefile.am
  echo *.cpp | fmt | sed -e 's/$/ \\/' -e '$s/\\$//' >> Makefile.am
  echo >>  Makefile.am
done

cd "$trunc_dir"
files=$(find . -type f -name "*.cpp" | grep -v detail)
[ -z "$files" ] && die "No cpp files in libs found"

echo -e "AM_CPPFLAGS = -I\$(srcdir)/..\n" >>  Makefile.am
echo -e "lib_LTLIBRARIES = libvotca_boost.la\n" >> Makefile.am
echo -e "libvotca_boost_la_SOURCES = \\" >> Makefile.am
echo -e "$files" | sed 's#^\./##' | fmt | sed -e 's/$/ \\/' -e '$s/\\$//' >> Makefile.am
echo >> Makefile.am
cd ..
echo Done with libs

echo Doing headers
cd boost
trunc_dir="$PWD"
files=$(find . -type f -not -name "*.hpp" -a -not -name "Makefile*")
[ -n "$files" ] && die "Unknown files:\n$files"

for i in $(find . -type d); do
  cd "$trunc_dir"
  echo Entering $i
  cd $i
  [ "$(echo *)" = '*' ] && continue
  rm -f Makefile.am
  dirs=""
  for j in *; do
    [ -d "$j" ] && dirs="$dirs $j" 
  done
  [ -n "$dirs" ] && echo -e "SUBDIRS = $dirs\n" >> Makefile.am
  [ "$(echo *.hpp)" = '*.hpp' ] && continue
  echo -e "boostincldir = \$(includedir)/votca/boost${i##.}\n" >> Makefile.am
  echo -e "boostincl_HEADERS = \\" >> Makefile.am
  echo *.hpp | fmt | sed -e 's/$/ \\/' -e '$s/\\$//' >> Makefile.am
  echo >>  Makefile.am
done
cd ${trunc_dir}
cd ..
echo Done with headers

echo Updating configure.ac
sed -n '1,/###ADD/p' configure.ac.in > configure.ac 
for i in $(find . -name Makefile.am); do
  echo $i | sed 's#^\./\(.*\).am#AC_CONFIG_FILES([\1])#' >> configure.ac
done
sed -n '/###END/,$p' configure.ac.in >> configure.ac
