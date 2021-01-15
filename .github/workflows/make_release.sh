#! /bin/bash -e

url="https://github.com/votca/votca.git"
branch=stable
testing=no
verbose=
what=( tools csg csg-tutorials xtp xtp-tutorials )
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

is_part() { #checks if 1st argument is part of the set given by other arguments
  [[ -z $1 || -z $2 ]] && die "${FUNCNAME[0]}: Missing argument"
  [[ " ${@:2} " = *" $1 "* ]]
}
export -f is_part

show_help() {
  cat << eof
This is the script to make release tarballs for VOTCA
Usage: ${0##*/} [OPTIONS] rel_version path/to/votca/checkout"
OPTIONS:
    --help          Show this help
    --test          Just test, do not commit stuff
    --branch BRANCH Use BRANCH instead of '$branch'
    --repos REPOS   Use repos instead of '${what[@]}'
    -j  JOBS        Jobs to use instead of '$j'
    --verbose       Do a verbose build
-D*                 Extra option to give to cmake 

Examples:  ${0##*/} --help
           ${0##*/} --test 1.2.3 srcdir

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
   --repos)
     what=( $2 )
     shift 2;;
   --branch)
     branch="$2"
     shift 2;;
   --test)
     testing=yes
     shift 1;;
   --verbose)
     verbose=yes
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

for i in tools csg csg-tutorials; do
  if ! is_part "$i" "${what[@]}"; then
    die "$i needs to be part of the repo selection"
  fi
done

[[ -z $2 ]] && die "${0##*/}: missing argument - no srcdir!\nTry ${0##*/} --help"
if [[ ${CI} != "true" && ${testing} = "no" && ${branch} != "stable" ]]; then
  die "branch ${branch} cannot be use without testing"
fi

shopt -s extglob

topdir="${PWD}"

rel="$1"
[[ ${CI} != "true" && $testing = "no" && ${rel} != 20???(.[1-9]|-rc.[1-9]) ]] && die "release has the wrong form"
srcdir="$2"
[[ -d $srcdir ]] || git clone --recursive "$url" "$srcdir"
pushd "${srcdir}"
srcdir="${PWD}"
[[ -f tools/CMakeLists.txt ]] || die "Checkout in $srcdir has no tools/CMakeLists.txt"
popd

instdir="${topdir}/install"
build="${topdir}/build"

set -e
cleanup() {
  [[ $testing = "no" ]] || return
  echo "####### ERROR ABOVE #########"
  pushd "${srcdir}"
  for p in . "${what[@]}"; do
    echo "$p"
    git -C "${p}" reset --hard "origin/${branch}" || true
    git -C "${p}" tag --delete "v${rel}" || true
  done
  popd
}
trap cleanup EXIT

pushd "${srcdir}"
git remote update --prune
git checkout "$branch" || die "Could not checkout $branch"
git pull --ff-only
for p in "${what[@]}"; do
  pushd "${p}"
  [[ -z "$(git ls-files -mo --exclude-standard)" ]] || die "There are modified or unknown files in $p"
  git remote update --prune
  git checkout "$branch" || die "Could not checkout $branch"
  git pull --ff-only
  [[ -z "$(git ls-files -mo --exclude-standard)" ]] || die "There are modified or unknown files in $p"
  if [[ $testing = "yes" ]]; then
    :
  elif [[ -f CMakeLists.txt ]]; then
    sed -i "/set(PROJECT_VERSION/s/\"[^\"]*\"/\"$rel\"/" CMakeLists.txt || die "sed of CMakeLists.txt failed"
    git add CMakeLists.txt
    if [[ -f CHANGELOG.rst ]]; then
      sed -i "/^Version ${rel}\>/s/released ..\...\.../released $(date +%d.%m.%y)/" CHANGELOG.rst
      git add CHANGELOG.rst
    fi
  fi
  if [[ $testing = "no" ]]; then
    [[ -f CHANGELOG.rst && -z $(grep "^Version ${rel}\>" CHANGELOG.rst) ]] && \
          die "Go and update CHANGELOG.rst in ${p} before making a release"
    #|| true because maybe version has not changed
    git commit -m "Version bumped to $rel" || true
    git tag "v${rel}"
  fi
  git archive --prefix "votca-${p}-${rel}/" -o "${topdir}/votca-${p}-${rel}.tar.gz" HEAD || die "git archive failed"
  popd
done
popd

rm -rf "$instdir" "$build"
mkdir "$instdir"
mkdir "$build"
pushd "$build"

echo "Starting build check from tarball"

cmake -DCMAKE_INSTALL_PREFIX="${instdir}" -DMODULE_BUILD=ON \
      -DVOTCA_TARBALL_DIR="${topdir}" -DVOTCA_TARBALL_TAG="${rel}" \
      -DENABLE_TESTING=ON \
      -DENABLE_REGRESSION_TESTING=ON \
      $(is_part xtp "${what[@]}" && echo -DBUILD_XTP=ON) \
      "${cmake_opts[@]}" "${srcdir}"
make -j"${j}" ${verbose:+VERBOSE=1}
popd

rm -rf "$build"
rm -rf "$instdir"

pushd "$srcdir"
if [[ $testing = "no" ]]; then
  sed -i "/set(PROJECT_VERSION/s/\"[^\"]*\"/\"$rel\"/" CMakeLists.txt || die "sed of CMakeLists.txt failed"
  if [[ -f README.rst ]]; then
    sed -i "/stable/s/or 'stable' or '[^']*'/or 'stable' or 'v$rel'/" README.rst share/doc/INSTALL.rst || die "sed of README.rst failed"
  fi
  git add -u
  git commit -m "Version bumped to $rel"
  git tag "v${rel}"
fi
trap - EXIT

if [[ $testing = "no" ]]; then
  echo "####### TODO by you #########"
  echo "cd $srcdir"
  echo "for p in . ${what[@]}; do git -C \$p log -p --submodule origin/${branch}..${branch}; done"
  echo "for p in . ${what[@]}; do git -C \$p  push --tags origin ${branch}:${branch}; done"
else
  echo "cd $topdir"
  echo "Take a look at " ./*"${rel}"*
fi
