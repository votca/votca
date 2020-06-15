#!/bin/bash -le

shopt -s extglob

die () {
  [[ -n $1 ]] && echo "$*" >&2
  exit 1
}

print_output() {
  [[ -n $1 ]] || die "${FUNCNAME[0]}: missing argument"
  echo "name=$1::${@:2}"
  echo "::set-output name=$1::${@:2}"
}

for i in INPUT_MINIMAL INPUT_OWN_GMX INPUT_REGRESSION_TESTING; do
  [[ "${!i}" = @(true|false) ]] || die "value of $i is ${!i}, excepted 'true' or 'false'"
  echo "$i='${!i}'"
done
[[ -n ${INPUT_DISTRO} ]] || die "value of INPUT_DISTRO was empty"
for i in INPUT_DISTRO INPUT_CMAKE_BUILD_TYPE INPUT_TOOLCHAIN INPUT_COVERAGE INPUT_CTEST_ARGS INPUT_CMAKE_ARGS; do
  echo "$i='${!i}'"
done

cmake_args=( -DCMAKE_VERBOSE_MAKEFILE=ON -DENABLE_TESTING=ON )
if [[ ${INPUT_CMAKE_BUILD_TYPE} ]]; then
  cmake_args+=( -DCMAKE_BUILD_TYPE=${INPUT_CMAKE_BUILD_TYPE} )
fi
if [[ ${INPUT_TOOLCHAIN} = "gnu" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc )
elif [[ ${INPUT_TOOLCHAIN} = "clang" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang -GNinja )
else
  die "Unknown INPUT_TOOLCHAIN"
fi
if [[ ${INPUT_COVERAGE} && ${INPUT_COVERAGE} != "false" ]]; then
  cmake_args+=( -DENABLE_COVERAGE_BUILD=ON )
  cov_tag=true
else
  cov_tag=false
fi

if ${INPUT_MODULE}; then
  cmake_args+=( -DMODULE_BUILD=ON -DCMAKE_INSTALL_PREFIX=$HOME/votca.install )
else
  cmake_args+=( -DCMAKE_INSTALL_PREFIX=/usr )
fi	
if ${INPUT_OWN_GMX}; then
  cmake_args+=( -DBUILD_OWN_GROMACS=ON -DENABLE_WARNING_FLAGS=OFF -DENABLE_WERROR=OFF )
else
  cmake_args+=( -DENABLE_WERROR=ON )
fi
if ${INPUT_MINIMAL}; then
  cmake_args+=( -DCMAKE_DISABLE_FIND_PACKAGE_HDF5=ON -DCMAKE_DISABLE_FIND_PACKAGE_FFTW3=ON -DCMAKE_DISABLE_FIND_PACKAGE_MKL=ON -DBUILD_MANPAGES=OFF -DCMAKE_DISABLE_FIND_PACKAGE_GROMACS=ON -DBUILD_XTP=OFF )
else
  cmake_args+=( -DBUILD_CSGAPPS=ON -DBUILD_XTP=ON -DBUILD_CSG_MANUAL=ON )
fi

if ${INPUT_MINIMAL} || [[ ${INPUT_DISTRO} = ubuntu@(|_rolling|_devel) ]];  then
  # Ubuntu 20.04 and above come with gromacs-2020, which doesn't have tabulated interaciton that are needed for csg regression tests
  # see https://gitlab.com/gromacs/gromacs/-/issues/1347
  # hopefully we can reenable this in the future with gromacs-2021
  cmake_args+=( -DENABLE_REGRESSION_TESTING=OFF )
else
  cmake_args+=( -DENABLE_REGRESSION_TESTING=${INPUT_REGRESSION_TESTING} )
fi

# lmp currently ill instruction on opensuse https://github.com/votca/csg-tutorials/issues/89
if [[ ${INPUT_DISTRO} = opensuse ]]; then
  cmake_args+=( -DCMAKE_DISABLE_FIND_PACKAGE_LMP=ON )
fi

cmake_args+=( ${INPUT_CMAKE_ARGS} )
print_output "cmake_args" "${cmake_args[@]}"

cache_key="ccache-${INPUT_DISTRO}-${INPUT_TOOLCHAIN}-${INPUT_CMAKE_BUILD_TYPE}-minimal-${INPUT_MINIMAL}-owngmx-${owngmx}-module-${INPUT_MODULE}-coverage-${cov_tag}"
print_output "cache_restore_key" "${cache_key}"
print_output "cache_key" "${cache_key}-$(date +%s)"

if [[ ${INPUT_DISTRO} = ubuntu*  ]] || ${INPUT_MODULE}; then
  # https://github.com/votca/votca/issues/318, sphinx build is currently broken on Ubuntu, due to sphinx 1.*
  # fedora uses sphinx 2.*
  print_output "build_sphinx" "false"
else
  print_output "build_sphinx" "true"
fi

if [[ ${INPUT_DISTRO} = "latest" ]] && ! ${INPUT_MODULE}; then
  print_output "check_format" "true"
else
  print_output "check_format" "false"
fi

# Grep project name from CMakeLists.txt and cut votca- suffix
[[ -f CMakeLists.txt ]] || die "No CMakeLists.txt found"
project=$(sed -n 's/project(\(votca-\)\?\([^)]*\))/\2/p' CMakeLists.txt)
[[ ${project} ]] || die "Could not fetch project"

ctest_args=( -L ${project} )
if [[ ${INPUT_COVERAGE} ]]; then
  # split coverage into 4 group with less than 1hr runtime
  # used votca/votca, csg, tools only
  # other modules can use 'RestGroup' to run all tests

  # false means the same as empty
  if [[ ${INPUT_COVERAGE} = "false" ]]; then
    :
  elif [[ ${INPUT_COVERAGE} = "Group1" ]]; then
    ctest_args+=( -R "regression_urea-water" )
  elif [[ ${INPUT_COVERAGE} = "Group2" ]]; then
    ctest_args+=( -R "'regression_spce_(re|imc|cma)'" )
  elif [[ ${INPUT_COVERAGE} = "Group3" ]]; then
    ctest_args+=( -R "'regression_(methanol-water|propane_imc)'" )
  elif [[ ${INPUT_COVERAGE} = "RestGroup" ]]; then
    ctest_args+=( -E "'regression_(urea-water|spce_(re|imc|cma)|methanol-water|propane_imc)'" )
  else
    die "Unknown coverage set"
  fi
fi
ctest_args+=( ${INPUT_CTEST_ARGS} )
print_output "ctest_args" "${ctest_args[@]}"

j="$(grep -c processor /proc/cpuinfo 2>/dev/null)" || j=0
((j++))
print_output "jobs" "${j}"

# Checkout votca main repo if we are building a module
if [[ ${project} != votca ]]; then
  git clone https://github.com/votca/votca
  if [[ ${GITHUB_REF} = refs/pull/*/merge ]]; then # pull request
    branch="${GITHUB_BASE_REF}"
  elif [[ ${GITHUB_REF} = refs/heads/* ]]; then # branch, e.g. stable
    branch=${GITHUB_REF#refs/heads/}
  elif [[ ${GITHUB_REF} = refs/tags/* ]]; then # tag or release
    branch=${GITHUB_REF#refs/tags/}
  else
    die "Handling on GITHUB_REF=${GITHUB_REF} not implemented"
  fi
  if [[ ${branch} && ${branch} != master ]]; then
    git -C votca checkout -b "${branch}" || true # || true as the branch might not exist
  fi
  git -C votca submodule update --init
  rm -rf "votca/${project}"
  ln -s "../../${project}" "votca/${project}"
fi
