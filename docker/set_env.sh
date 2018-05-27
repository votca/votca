#!/bin/bash -xe

die () {
  echo "$*" >&2
  exit 1
}

set -x
if [[ $ENV -eq 1 ]]; then
  # Debug build with half the tests (due to timing constraints)
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)"
  export CMAKE_BUILD_TYPE=Debug
elif [[ $ENV -eq 2 ]]; then
  # Debug build with second half of the tests
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -R _re"
  export CMAKE_BUILD_TYPE=Debug
elif [[ $ENV -eq 3 ]]; then
  # Debug build with -Werror, build tests as well, but exclude all tests
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E ."
  export CMAKE_BUILD_TYPE=Debug
  export WERROR=yes
elif [[ $ENV -eq 4 ]]; then
  # Release build, which gets push to dockerhun
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  export CMAKE_BUILD_TYPE=Release
  export DOCKERHUB=yes
elif [[ $ENV -eq 5 ]]; then
  # Release build with gromacs-2016
  export DISTRO=fedora_gmx2016
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  export CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 6 ]]; then
  # Release build with gromacs-2016 (double)
  export DISTRO=fedora_gmx2016_d
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  export CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 7 ]]; then
  # Release build with gromacs-2018
  export DISTRO=fedora_gmx2018
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  export CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 8 ]]; then
  # Release build with gromacs-2018 (double)
  export DISTRO=fedora_gmx2018_d
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  export CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 9 ]]; then
  # Release build with gromacs master
  export DISTRO=fedora_gmx9999
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  export CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 10 ]]; then
  # Release build with gromacs master (double)
  export DISTRO=fedora_gmx9999_d
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  export CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 11 ]]; then
  # Release build on Ubuntu
  export DISTRO=ubuntu
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  export CMAKE_BUILD_TYPE=Release
elif [[ $ENV -eq 12 ]]; then
  # Release build with -Werror, build tests as well, but exclude all tests
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E ."
  export CMAKE_BUILD_TYPE=Release
  export WERROR=yes
elif [[ $ENV -eq 13 ]]; then
  # Build with no cmake_build_type and coverage on, first half of the tests
  # superbuild has no code, so no coverage, but test None build type
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)"
  export CMAKE_BUILD_TYPE=None
  export COVERAGE=yes
  export SKIP=yes # bug #67
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
elif [[ $ENV -eq 14 ]]; then
  # Build with no cmake_build_type and coverage on, second half of the tests
  # superbuild has no code, so no coverage, but test None build type
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -R _re"
  export CMAKE_BUILD_TYPE=None
  export COVERAGE=yes
  export SKIP=yes # bug #67
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
elif [[ $ENV -eq 15 ]]; then
  # Build with doxygen
  export TESTING=OFF
  export CMAKE_BUILD_TYPE=Release
  export DOXYGEN=yes
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
  [[ ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes #only votca/votca actually deploys doxygen - re-enable if we check for doxygen warnings
elif [[ $ENV -eq 16 ]]; then
  # minimal build without external libs (tools & csg only)
  export TESTING=OFF
  export CMAKE_BUILD_TYPE=Release
  export MINIMAL=yes
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */tools || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # minimal build is only for tools, csg and superbuild
elif [[ $ENV -eq 17 ]]; then
  # module build
  export MODULE_BUILD=ON
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)"
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
elif [[ $ENV -eq 18 ]]; then
  # build internal gromacs
  export DISTRO=fedora_nogmx
  export BUILD_GROMACS=ON
  export TESTING=ON
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)"
  export CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
else
  die "Unknown environment"
fi
set +x
