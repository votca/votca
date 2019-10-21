#!/bin/bash -xe

shopt -s extglob

die () {
  echo "$*" >&2
  exit 1
}


docker_opts=()
add_to_docker_opts() {
  docker_opts+=( "--build-arg" "$1" )
}

set -x

for i in CC CXX CI TRAVIS TRAVIS_BRANCH TRAVIS_JOB_NUMBER TRAVIS_PULL_REQUEST TRAVIS_JOB_ID TRAVIS_TAG TRAVIS_REPO_SLUG \
  TRAVIS_COMMIT TRAVIS_PULL_REQUEST_SHA; do
  [[ -z ${!i} ]] || add_to_docker_opts "$i=${!i}"
done

if [[ $ENV -eq 1 ]]; then
  # Release build with doxygen, which gets push to dockerhub, first half of the tests
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E _re -LE memory"
  # only run csg-tutorials regressions tests for csg, tools, csg-tutorials and votca
  [[ ${TRAVIS_REPO_SLUG#*/} = @(csg|tools|csg-tutorials|votca) && ${TRAVIS_BUILD_STAGE_NAME} != "Deploy" ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  export WERROR=yes
  add_to_docker_opts DOXYGEN=yes
  [[ ${TRAVIS_BUILD_STAGE_NAME} != "Deploy" ]] && add_to_docker_opts DOXYGEN_COVERAGE=yes
elif [[ $ENV -eq 2 ]]; then
  # Release build, which gets push to dockerhub, first half of the tests (for tools, csg, csg-tutorials and votca)
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -R _re -LE memory"
  # only run csg-tutorials regressions tests for csg, tools, csg-tutorials and votca
  [[ ${TRAVIS_REPO_SLUG#*/} = @(csg|tools|csg-tutorials|votca) ]] && add_to_docker_opts REGRESSION_TESTING=ON || export SKIP=yes
  add_to_docker_opts CMAKE_BUILD_TYPE=Release
  export WERROR=yes
elif [[ $ENV -eq 3 ]]; then
  # Debug build (no time problem as no tests are run)
  add_to_docker_opts TESTING=OFF
  add_to_docker_opts CMAKE_BUILD_TYPE=Debug
  add_to_docker_opts CLANG_FORMAT=yes
  export WERROR=yes
elif [[ $ENV -eq 4 ]]; then
  # Coverage build on Ubuntu, Fedora has issues, see #67
  # Build with no cmake_build_type and coverage on, first half of the tests
  # superbuild has no code, so no coverage, but test "None" build type
  export DISTRO=ubuntu
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E _re -LE memory"
  # only run csg-tutorials regressions tests for csg, tools - csg-tutorials and votca have no coverable code
  [[ ${TRAVIS_REPO_SLUG#*/} = @(csg|tools) ]] && add_to_docker_opts REGRESSION_TESTING=ON
  add_to_docker_opts CMAKE_BUILD_TYPE=
  CXXFLAGS="-O2" # gcc's default would be -O0 (slow!) otherwise
  add_to_docker_opts COVERAGE=yes
  [[ $CXX = clang++ ]] && export SKIP=yes #clang coverage is too slow for travis, plus we have gcc's coverage, so not much more insights
elif [[ $ENV -eq 5 ]]; then
  # Coverage build on Ubuntu, Fedora has issues, see #67
  # Build with no cmake_build_type and coverage on, second half of the tests (for tools and csg)
  # superbuild has no code, so no coverage, but test "None" build type
  export DISTRO=ubuntu
  add_to_docker_opts TESTING=ON
  add_to_docker_opts TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -R _re -LE memory"
  # only run csg-tutorials regressions tests for csg, tools - csg-tutorials and votca have no coverable code
  [[ ${TRAVIS_REPO_SLUG#*/} = @(csg|tools) ]] && add_to_docker_opts REGRESSION_TESTING=ON || export SKIP=yes
  add_to_docker_opts CMAKE_BUILD_TYPE=
  CXXFLAGS="-O2" # gcc's default would be -O0 (slow!) otherwise
  add_to_docker_opts COVERAGE=yes
  [[ $CXX = clang++ ]] && export SKIP=yes #clang coverage is too slow for travis, plus we have gcc's coverage, so not much more insights
elif [[ $ENV -eq 6 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 7 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 8 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 9 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 10 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 11 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 12 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 13 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 14 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 15 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 16 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 17 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 18 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 19 ]]; then
  export SKIP=yes
elif [[ $ENV -eq 20 ]]; then
  export SKIP=yes
else
  die "Unknown environment"
fi

add_to_docker_opts CXXFLAGS="${CXXFLAGS} -Wall -Wextra -Wpedantic -Wshadow -Wconversion ${WERROR:+-Werror}"
add_to_docker_opts TRAVIS_OS_NAME="${DISTRO:-fedora}"

export docker_opts
set +x
