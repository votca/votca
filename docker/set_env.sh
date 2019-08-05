#!/bin/bash -xe

die () {
  echo "$*" >&2
  exit 1
}


docker_opts=()
add_to_docker_opts() {
  docker_opts+=( "--build-arg" "$1" )
}

set -x
if [[ $ENV -eq 1 ]]; then
  # Debug build with half the tests (due to timing constraints)
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Debug
elif [[ $ENV -eq 2 ]]; then
  # Debug build with second half of the tests
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -R _re"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Debug
elif [[ $ENV -eq 3 ]]; then
  # Debug build with -Werror, build tests as well, but exclude all tests
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E ."
  add_to_docker_opts CMAKE_BUILD_TYPE=Debug
  export WERROR=yes
elif [[ $ENV -eq 4 ]]; then
  # Release build, which gets push to dockerhub
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  export DOCKERHUB=yes
elif [[ $ENV -eq 5 ]]; then
  # Release build with gromacs-2016
  export DISTRO=fedora_gmx2016
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 6 ]]; then
  # Release build with gromacs-2016 (double)
  export DISTRO=fedora_gmx2016_d
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 7 ]]; then
  # Release build with gromacs-2018
  export DISTRO=fedora_gmx2018
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 8 ]]; then
  # Release build with gromacs-2018 (double)
  export DISTRO=fedora_gmx2018_d
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 9 ]]; then
  # Release build with gromacs master
  export DISTRO=fedora_gmx9999
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
  export SKIP=yes # bug votca/csg#387
elif [[ $ENV -eq 10 ]]; then
  # Release build with gromacs master (double)
  export DISTRO=fedora_gmx9999_d
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
  export SKIP=yes # bug votca/csg#387
elif [[ $ENV -eq 11 ]]; then
  # Release build on Ubuntu
  export DISTRO=ubuntu
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
elif [[ $ENV -eq 12 ]]; then
  # Release build with -Werror, build tests as well, but exclude all tests
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E ."
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  export WERROR=yes
elif [[ $ENV -eq 13 ]]; then
  # Build with no cmake_build_type and coverage on, first half of the tests
  # superbuild has no code, so no coverage, but test None build type
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=None
  add_to_docker_opts COVERAGE=yes
  export SKIP=yes # bug #67
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
elif [[ $ENV -eq 14 ]]; then
  # Build with no cmake_build_type and coverage on, second half of the tests
  # superbuild has no code, so no coverage, but test None build type
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -R _re"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=None
  add_to_docker_opts COVERAGE=yes
  export SKIP=yes # bug #67
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
elif [[ $ENV -eq 15 ]]; then
  # Build with doxygen
  add_to_docker_opts TESTING=OFF
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  add_to_docker_opts DOXYGEN=yes
  export DOXYGEN=yes # uses .travis.yml
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
  [[ ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes #only votca/votca actually deploys doxygen - re-enable if we check for doxygen warnings
elif [[ $ENV -eq 16 ]]; then
  # minimal build without external libs (tools & csg only)
  export DISTRO=fedora_nogmx
  add_to_docker_opts TESTING=OFF
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  add_to_docker_opts MINIMAL=yes
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */tools || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # minimal build is only for tools, csg and superbuild
elif [[ $ENV -eq 17 ]]; then
  # module build
  add_to_docker_opts MODULE_BUILD=ON
  [[ $CC = clang ]] && add_to_docker_opts TESTING=ON || add_to_docker_opts TESTING=OFF
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
elif [[ $ENV -eq 18 ]]; then
  # build internal gromacs
  export DISTRO=fedora_nogmx
  add_to_docker_opts BUILD_GROMACS=ON
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
elif [[ $ENV -eq 19 ]]; then
  export DISTRO=fedora_gmx2019
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
elif [[ $ENV -eq 20 ]]; then
  export DISTRO=fedora_gmx2019_d
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)"
  [[ ${TRAVIS_REPO_SLUG} = */csg-tutorials ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  [[ ${TRAVIS_REPO_SLUG} = */csg || ${TRAVIS_REPO_SLUG} = */votca ]] || export SKIP=yes # only csg uses gromacs
else
  die "Unknown environment"
fi
export docker_opts
set +x
