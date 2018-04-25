#!/bin/bash -xe

die () {
  echo "$*" >&2
  exit 1
}

if [[ $ENV -eq 1 ]]; then 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)" 
  export CMAKE_BUILD_TYPE=Debug
elif [[ $ENV -eq 2 ]]; then 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -R _re" 
  export CMAKE_BUILD_TYPE=Debug
elif [[ $ENV -eq 3 ]]; then 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E ." 
  export CMAKE_BUILD_TYPE=Debug   
  export WERROR=yes
elif [[ $ENV -eq 4 ]]; then 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)" 
  export CMAKE_BUILD_TYPE=Release 
  export DOCKERHUB=yes
elif [[ $ENV -eq 5 ]]; then 
  export DISTRO=fedora_gmx2016 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)" 
  export CMAKE_BUILD_TYPE=Release
elif [[ $ENV -eq 6 ]]; then 
  export DISTRO=fedora_gmx2016_d 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)" 
  export CMAKE_BUILD_TYPE=Release
elif [[ $ENV -eq 7 ]]; then 
  export DISTRO=fedora_gmx2018 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)" 
  export CMAKE_BUILD_TYPE=Release
elif [[ $ENV -eq 8 ]]; then 
  export DISTRO=fedora_gmx2018_d 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)" 
  export CMAKE_BUILD_TYPE=Release
elif [[ $ENV -eq 9 ]]; then 
  export DISTRO=fedora_gmx9999 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)" 
  export CMAKE_BUILD_TYPE=Release
elif [[ $ENV -eq 10 ]]; then 
  export DISTRO=fedora_gmx9999_d 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)" 
  export CMAKE_BUILD_TYPE=Release
elif [[ $ENV -eq 11 ]]; then 
  export DISTRO=ubuntu 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\)" 
  export CMAKE_BUILD_TYPE=Release
elif [[ $ENV -eq 12 ]]; then 
  export TESTING=OFF 
  export CMAKE_BUILD_TYPE=Release 
  export WERROR=yes
elif [[ $ENV -eq 13 ]]; then 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -E \(_imc\|spce_cma_simple\|_re\)" 
  export CMAKE_BUILD_TYPE=None 
#useless for votca/votca, but for coverage of individual modules
  export COVERAGE=yes
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
elif [[ $ENV -eq 14 ]]; then 
  export TESTING=ON 
  export TESTOPTS="-L ${TRAVIS_REPO_SLUG#*/} -R _re" 
  export CMAKE_BUILD_TYPE=None 
  export COVERAGE=yes
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
elif [[ $ENV -eq 15 ]]; then 
  export TESTING=OFF 
  export CMAKE_BUILD_TYPE=Release 
  export DOXYGEN=yes
  [[ $CC = clang ]] && export SKIP=yes # no new info when using clang
elif [[ $ENV -eq 16 ]]; then 
  export TESTING=OFF 
  export CMAKE_BUILD_TYPE=Release 
  export MINIMAL=yes
else
  die "Unknown enviorment"
fi
