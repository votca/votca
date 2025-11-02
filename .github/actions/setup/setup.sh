#!/bin/bash -le

shopt -s extglob

die () {
  [[ -n $1 ]] && echo "$*" >&2
  exit 1
}

print_output() {
  [[ -n $1 ]] || die "${FUNCNAME[0]}: missing argument"
  echo "Setting $1=${@:2}"
  echo "$1=${@:2}" >> $GITHUB_OUTPUT
}

for i in INPUT_MINIMAL; do
  [[ "${!i}" = @(true|false) ]] || die "value of $i is ${!i}, excepted 'true' or 'false'"
  echo "$i='${!i}'"
done
[[ -n ${INPUT_DISTRO} ]] || die "value of INPUT_DISTRO was empty"
for i in INPUT_DISTRO INPUT_CMAKE_BUILD_TYPE INPUT_TOOLCHAIN INPUT_COVERAGE INPUT_CTEST_ARGS INPUT_CMAKE_ARGS INPUT_CODE_ANALYZER; do
  echo "$i='${!i}'"
done

if [[ ${INPUT_BRANCH} ]]; then # user overwrite
  branch="${INPUT_BRANCH}"
elif [[ ${GITHUB_REF} = refs/pull/*/merge ]]; then # pull request
  branch="${GITHUB_BASE_REF}"
elif [[ ${GITHUB_REF} = refs/heads/* ]]; then # branch, e.g. stable
  branch=${GITHUB_REF#refs/heads/}
elif [[ ${GITHUB_REF} = refs/tags/* ]]; then # tag or release
  branch=${GITHUB_REF#refs/tags/}
else
  die "Handling on GITHUB_REF=${GITHUB_REF} not implemented"
fi

cmake_args=( -DCMAKE_VERBOSE_MAKEFILE=ON -DINSTALL_CSGAPPS=ON )
# do not inject -march=native as the CI runs on different backends and hence will create a conflict with ccache
cmake_args+=( -DINJECT_MARCH_NATIVE=OFF )
if [[ ${INPUT_CMAKE_BUILD_TYPE} ]]; then
  cmake_args+=( -DCMAKE_BUILD_TYPE=${INPUT_CMAKE_BUILD_TYPE} )
fi
if [[ ${INPUT_TOOLCHAIN} = "gnu" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=g++ )
elif [[ ${INPUT_TOOLCHAIN} = "clang" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=clang++ -GNinja )
elif [[ ${INPUT_TOOLCHAIN} = "intel" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=icpc )
elif [[ ${INPUT_TOOLCHAIN} = "intel-oneapi" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=icpx )
elif [[ ${INPUT_TOOLCHAIN} = "intel-oneapi-dpc" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=dpcpp )
else
  die "Unknown INPUT_TOOLCHAIN; ${INPUT_TOOLCHAIN}"
fi

cmake_args+=( -DCMAKE_INSTALL_PREFIX=/usr )

# use ccache
cmake_args+=( -DCMAKE_CXX_COMPILER_LAUNCHER=ccache )

cmake_args+=( -DENABLE_WERROR=ON )
if [[ ${INPUT_MINIMAL} = true ]]; then
  cmake_args+=( -DCMAKE_DISABLE_FIND_PACKAGE_HDF5=ON -DCMAKE_DISABLE_FIND_PACKAGE_FFTW3=ON -DCMAKE_DISABLE_FIND_PACKAGE_MKL=ON -DCMAKE_DISABLE_FIND_PACKAGE_GROMACS=ON -DBUILD_MANPAGES=OFF -DBUILD_XTP=OFF -DBUILD_TESTING=OFF )
else
  cmake_args+=( -DBUILD_XTP=ON )
fi

if [[ ${INPUT_DISTRO} = *:intel ]]; then
  cmake_args+=( -DREQUIRE_MKL=ON )
fi

if [[ ${INPUT_CODE_ANALYZER} = "codeql" ]]; then
  # CodeQL does not work with valgrind
  cmake_args+=( -DVALGRIND_EXECUTABLE=FALSE )
elif [[ ${INPUT_COVERAGE} && ${INPUT_COVERAGE} != "false" ]] || [[ ${INPUT_CODE_ANALYZER} = coverage* ]]; then
  cmake_args+=( -DENABLE_COVERAGE_BUILD=ON )
elif [[ ${INPUT_CODE_ANALYZER} = "false" ]]; then
  :
else
  die "Unknown INPUT_CODE_ANALYZER: ${INPUT_CODE_ANALYZER}"
fi

cmake_args+=( ${INPUT_CMAKE_ARGS} )
print_output "cmake_args" "${cmake_args[@]}"

cache_key="ccache-${INPUT_DISTRO/:/_}-${INPUT_TOOLCHAIN}-${INPUT_CMAKE_BUILD_TYPE}-minimal-${INPUT_MINIMAL}-analysis-${INPUT_CODE_ANALYZER%%:*}"
print_output "cache_restore_key" "${cache_key}"
print_output "cache_key" "${cache_key}-$(date +%s)"

if [[ ${branch} = stable || ${INPUT_DISTRO} != fedora:@(latest|rawhide)  || ${INPUT_CMAKE_BUILD_TYPE} = Debug ]]; then
  # Only build doc sphinx on Fedora, as there many issues on other, e.g.:
  # 1.) Don't build sphinx on stable, not useful, only master is useful
  # 2.) On Ubuntu 18.04 sphinx is too old for nbsphinx
  #     File "/usr/lib/python3/dist-packages/nbsphinx.py", line 1383, in _add_notebook_parser
  #       source_suffix.append('.ipynb')
  #     AttributeError: 'dict' object has no attribute 'append'
  #     nbsphinx that requires that sphinx>1.8 but in Ubuntu 18.04 sphinx==1.6.7
  # 3.) Debug builds are too slow to run notebooks in xtp-tutorials
  print_output "build_sphinx" "false"
else
  print_output "build_sphinx" "true"
fi

if [[ ${INPUT_DISTRO} = "fedora:latest" ]]; then
  print_output "check_format" "true"
else
  print_output "check_format" "false"
fi

ctest_args=( )
ctest_args+=( ${INPUT_CTEST_ARGS} )
print_output "ctest_args" "${ctest_args[@]}"

j="$(grep -c processor /proc/cpuinfo 2>/dev/null)" || j=0
((j++))
print_output "jobs" "${j}"
