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
for i in INPUT_DISTRO INPUT_CMAKE_BUILD_TYPE INPUT_TOOLCHAIN INPUT_COVERAGE INPUT_CTEST_ARGS INPUT_CMAKE_ARGS INPUT_MODULE; do
  echo "$i='${!i}'"
done

# Grep module name from CMakeLists.txt and cut votca- suffix
[[ -f CMakeLists.txt ]] || die "No CMakeLists.txt found"
module=$(sed -n 's/project(\(votca-\)\?\([^ )]*\).*)/\2/p' CMakeLists.txt)
[[ ${module} ]] || die "Could not fetch module"
print_output "module" "${module}"

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

cmake_args=( -DCMAKE_VERBOSE_MAKEFILE=ON -DENABLE_TESTING=ON  -DBUILD_CSGAPPS=ON )
if [[ ${INPUT_CMAKE_BUILD_TYPE} ]]; then
  cmake_args+=( -DCMAKE_BUILD_TYPE=${INPUT_CMAKE_BUILD_TYPE} )
fi
if [[ ${INPUT_TOOLCHAIN} = "gnu" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=g++ )
elif [[ ${INPUT_TOOLCHAIN} = "clang" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=clang++ -GNinja )
elif [[ ${INPUT_TOOLCHAIN} = "intel" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=icpc )
  mkdir ~/Licenses
  curl https://dynamicinstaller.intel.com/api/v2/license > ~/Licenses/intel.lic
elif [[ ${INPUT_TOOLCHAIN} = "intel-oneapi" ]]; then
  cmake_args+=( -DCMAKE_CXX_COMPILER=icpx )
else
  die "Unknown INPUT_TOOLCHAIN"
fi

if [[ ${INPUT_COVERAGE} && ${INPUT_COVERAGE} != "false" ]]; then
  cmake_args+=( -DENABLE_COVERAGE_BUILD=ON )
  cov_tag=true
else
  cov_tag=false
fi

if [[ ${INPUT_MODULE} = true ]]; then
  cmake_args+=( -DMODULE_BUILD=ON -DCMAKE_INSTALL_PREFIX=$HOME/votca.install )
else
  cmake_args+=( -DCMAKE_INSTALL_PREFIX=/usr )
fi	
if [[ ${INPUT_OWN_GMX} = true ]]; then
  cmake_args+=( -DBUILD_OWN_GROMACS=ON -DENABLE_WARNING_FLAGS=OFF -DENABLE_WERROR=OFF )
  # remove this block when gromacs uses cxx only, i.e. gmx2021
  if [[ ${INPUT_TOOLCHAIN} = "gnu" ]]; then
    cmake_args+=( -DCMAKE_C_COMPILER=gcc )
  elif [[ ${INPUT_TOOLCHAIN} = "clang" ]]; then
    cmake_args+=( -DCMAKE_C_COMPILER=clang )
  elif [[ ${INPUT_TOOLCHAIN} = "intel" ]]; then
    cmake_args+=( -DCMAKE_C_COMPILER=icc )
  elif [[ ${INPUT_TOOLCHAIN} = "intel-oneapi" ]]; then
    cmake_args+=( -DCMAKE_CXX_COMPILER=icx )
  fi
else
  cmake_args+=( -DENABLE_WERROR=ON )
fi
if [[ ${INPUT_MINIMAL} = true ]]; then
  cmake_args+=( -DCMAKE_DISABLE_FIND_PACKAGE_HDF5=ON -DCMAKE_DISABLE_FIND_PACKAGE_FFTW3=ON -DCMAKE_DISABLE_FIND_PACKAGE_MKL=ON -DCMAKE_DISABLE_FIND_PACKAGE_GROMACS=ON -DBUILD_MANPAGES=OFF -DBUILD_XTP=OFF )
elif [[ ${module} = csg-tutorials ]]; then
  cmake_args+=( -DBUILD_XTP=OFF )
else
  cmake_args+=( -DBUILD_XTP=ON )
fi

if [[ ${INPUT_MINIMAL} = true || ${INPUT_DISTRO} = ubuntu:@(latest|rolling|devel) ]];  then
  # Ubuntu 20.04 and above come with gromacs-2020, which doesn't have tabulated interaciton that are needed for csg regression tests
  # see https://gitlab.com/gromacs/gromacs/-/issues/1347
  # hopefully we can reenable this in the future with gromacs-2021
  cmake_args+=( -DENABLE_REGRESSION_TESTING=OFF )
else
  cmake_args+=( -DENABLE_REGRESSION_TESTING=${INPUT_REGRESSION_TESTING} )
fi

cmake_args+=( ${INPUT_CMAKE_ARGS} )
print_output "cmake_args" "${cmake_args[@]}"

cache_key="ccache-${INPUT_DISTRO/:/_}-${INPUT_TOOLCHAIN}-${INPUT_CMAKE_BUILD_TYPE}-minimal-${INPUT_MINIMAL}-owngmx-${owngmx}-module-${INPUT_MODULE}-coverage-${cov_tag}"
print_output "cache_restore_key" "${cache_key}"
print_output "cache_key" "${cache_key}-$(date +%s)"

if [[ ${branch} = stable || ${INPUT_DISTRO} != fedora:@(latest|rawhide)  || ${INPUT_CMAKE_BUILD_TYPE} = Debug || ${INPUT_MODULE} = true ]]; then
  # Only build doc sphinx on Fedora, as there many issues on other, e.g.:
  # 1.) Don't build sphinx on stable, not useful, only master is useful
  # 2.) On Ubuntu 18.04 sphinx is too old for nbsphinx
  #     File "/usr/lib/python3/dist-packages/nbsphinx.py", line 1383, in _add_notebook_parser
  #       source_suffix.append('.ipynb')
  #     AttributeError: 'dict' object has no attribute 'append'
  #     nbsphinx that requires that sphinx>1.8 but in Ubuntu 18.04 sphinx==1.6.7
  # 3.) Debug builds are too slow to run notebooks in xtp-tutorials
  # 4.) Module build doesn't support sphinx
  print_output "build_sphinx" "false"
else
  print_output "build_sphinx" "true"
fi

if [[ ${INPUT_DISTRO} = "fedora:latest" ]] && [[ ${INPUT_MODULE} = false ]]; then
  print_output "check_format" "true"
else
  print_output "check_format" "false"
fi

ctest_args=( -L ${module} )
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
  elif [[ ${INPUT_COVERAGE} = "RestGroup" || ${INPUT_COVERAGE} = "true" ]]; then
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
if [[ ${module} != votca ]]; then
  git clone https://github.com/votca/votca
  if [[ ${branch} && ${branch} != master ]]; then
    git -C votca checkout "${branch}" || true # || true as the branch might not exist
  fi
  git -C votca submodule update --init
  git -C "votca/${module}" fetch "$PWD"
  git -C "votca/${module}" checkout FETCH_HEAD
fi
