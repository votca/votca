#! /bin/bash -e

url="https://github.com/votca/votca.git"
branch=master
verbose=
cmake_opts=()

die () {
  echo -e "$*"
  exit 1
}

unset CSGSHARE VOTCASHARE

# So that realtime test don't try to run X11
export GNUTERM=dumb

j="$(grep -c processor /proc/cpuinfo 2>/dev/null)" || j=0
((j++))

show_help() {
  cat << eof
This is the script to make release tarballs for VOTCA
Usage: ${0##*/} [OPTIONS] rel_version path/to/votca/checkout"
OPTIONS:
    --help          Show this help
    --branch BRANCH Use BRANCH instead of '$branch'
    -j  JOBS        Jobs to use instead of '$j'
    --verbose       Do a verbose build
    --debug         Run in debug more
-D*                 Extra option to give to cmake

Examples:  ${0##*/} --help
           ${0##*/} 1.2.3 srcdir

Report bugs and comments at https://github.com/votca/admin/issues
eof
}

shopt -s extglob
while [[ $# -gt 0 ]]; do
  if [[ ${1} = --*=* ]]; then # case --xx=yy
    set -- "${1%%=*}" "${1#*=}" "${@:2}" # --xx=yy to --xx yy
  elif [[ ${1} = -[^-]?* ]]; then # case -xy split
    if [[ ${1} = -[jpD]* ]]; then #short opts with arguments
       set -- "${1:0:2}" "${1:2}" "${@:2}" # -xy to -x y
    else #short opts without arguments
       set -- "${1:0:2}" "-${1:2}" "${@:2}" # -xy to -x -y
    fi
 fi
 case $1 in
   --branch)
     branch="$2"
     shift 2;;
   --verbose)
     verbose=yes
     shift 1;;
   --debug)
     set -x
     shift 1;;
   -D)
    cmake_opts+=( -D"${2}" )
    shift 2;;
   -j)
    j="$2"
    shift 2;;
   --help)
     show_help
     exit $?;;
   -*)
     die "Unknown options $1";;
   --)
     break;;
   *)
     break;;
 esac
done

[[ -z $2 ]] && die "${0##*/}: missing argument - no srcdir!\nTry ${0##*/} --help"

shopt -s extglob

topdir="${PWD}"

rel="$1"
[[ ${CI} != "true" && ${rel} != 20???(.[1-9]|-rc.[1-9]) ]] && die "release has the wrong form"
srcdir="$2"
[[ -d $srcdir ]] || git clone --recursive "$url" "$srcdir"
pushd "${srcdir}"
srcdir="${PWD}"
[[ -f tools/CMakeLists.txt ]] || die "Checkout in $srcdir has no tools/CMakeLists.txt"
popd

instdir="${topdir}/install"
build="${topdir}/build"

set -e

pushd "${srcdir}"
git remote update --prune
git checkout "$branch" || die "Could not checkout $branch"
git pull --ff-only
[[ -z "$(git ls-files -mo --exclude-standard)" ]] || die "There are modified or unknown files"

cur_rel="$(sed -n 's/set(PROJECT_VERSION *"\([^"]*\)").*/\1/p' CMakeLists.txt)"
[[ -n ${cur_rel} ]] || die "Could not grep current release from CMakeLists.txt)"

# check if CHANGELOG section has no entry, there should be at least something like "-  no changes"
version_section="$(awk -v r="^Version ${cur_rel}( |$)" '($0 ~ "^Version"){go=0} ($0 ~ r){go=1}{if(go==1){print $0}}' CHANGELOG.rst)"
line_nr="$(sed -n "/^Version ${cur_rel}\( \|$\)/=" CHANGELOG.rst)"
[[ $version_section ]] || die "Could not find section to $cur_rel"
echo "Found section for $cur_rel (starting line ${line_nr})"
last_line="$(echo "$version_section" | sed '/^[[:space:]]*$/d' | sed -n '$p')"
[[ $last_line ]] || die "Could not grep last line"
[[ ${last_line} = -* || ${last_line} = '   '[^\ ]* ]] || die "Last line isn't an item (does not start with -), but ${last_line}, fix the CHANGELOG.rst in $p first"

sed -i "/set(PROJECT_VERSION/s/\"[^\"]*\"/\"$rel\"/" CMakeLists.txt || die "sed of CMakeLists.txt failed"
git add CMakeLists.txt
new_section="Version ${rel} (released $(date +%d.%m.%y))"
sed -i "$((line_nr+1))d" CHANGELOG.rst || die "sed to delete === line failed"
sed -i "/^Version ${cur_rel}/s/.*/${new_section}\n${new_section//?/=}/" CHANGELOG.rst || die "sed of CHANGELOG.rst failed"
git add CHANGELOG.rst
sed -i "/version=/s/master # or '[^']*'/master # or 'v$rel'/" README.rst || die "sed of README.rst failed"
git add README.rst
sed -i "/version=/s/master # or '[^']*'/master # or 'v$rel'/" share/sphinx/INSTALL.rst || die "sed of INSTALL.rst failed"
git add share/sphinx/INSTALL.rst

#|| true because maybe version has not changed
git commit -m "Version bumped to $rel" || true
git tag "v${rel}"

git archive --prefix "votca-${rel}/" -o "${topdir}/votca-${rel}.tar.gz" HEAD || die "git archive failed"
popd

rm -rf "$instdir" "$build"
mkdir "$instdir"

echo "Starting build check from tarball"

tar -xvf "${topdir}/votca-${rel}.tar.gz"
cmake -DCMAKE_INSTALL_PREFIX="${instdir}" \
      -DBUILD_XTP=ON \
      -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
      "${cmake_opts[@]}" -S "votca-${rel}/" -B "$build"
cmake --build "${build}" -j"${j}" ${verbose:+--verbose}

rm -rf "$instdir" "$build"

pushd "${srcdir}"
add_rel="${rel}-dev"

sed -i "/set(PROJECT_VERSION/s/\"[^\"]*\"/\"${add_rel}\"/" CMakeLists.txt || die "sed of CMakeLists.txt failed"
git add CMakeLists.txt
git commit -m "Version bumped to ${add_rel}" || true

new_section="Version ${add_rel}"
sed -i "/^Version ${rel}\( \|$\)/i ${new_section}\n${new_section//?/=}\n" CHANGELOG.rst
git add CHANGELOG.rst
git commit -m "CHANGELOG: add ${add_rel} section"
popd

echo "####### TODO by you #########"
echo "cd $srcdir"
echo "git -C \$p log -p --submodule origin/${branch}..${branch}"
echo "git -C \$p  push --tags origin ${branch}:${branch}"
