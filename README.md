[![Codacy Badge](https://app.codacy.com/project/badge/Grade/b5567bfcf2c8411a8057c47fa7126781)](https://www.codacy.com/gh/votca/votca?utm_source=github.com&utm_medium=referral&utm_content=votca/votca&utm_campaign=Badge_Grade)
[![CI](https://github.com/votca/votca/workflows/CI/badge.svg?branch=stable)](https://github.com/votca/votca/actions?query=workflow%3ACI+branch%3Astable)
[![Docker](https://github.com/votca/votca/workflows/Docker/badge.svg?branch=stable)](https://github.com/votca/votca/actions?query=workflow%3A+branch%3Astable)

This is VOTCA's next generation build system for CSG and XTP. It allows you to easily install: 

-   VOTCA-CSG, a library which provides tools to develop coarse-grained potentials from atomistic simulation data
-   VOTCA-XTP, a library designed to determine electronic properties of organic materials from atomistic MD-trajectories.

Usage:

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'stable' or 'v1.6.2'
    git clone -b ${version} --recursive https://github.com/votca/votca.git
    cd votca
    mkdir build
    cd build
    cmake -DBUILD_CSGAPPS=ON -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} ..
    cmake --build . -- -j<number of cores>
    cmake --build . --target install

Using this code via docker:

    docker pull votca/votca
    docker run -it votca/votca /bin/bash

For further details see:

1. [Installation](share/doc/INSTALL.md)
2. [Further Information](http://www.votca.org)
3. [Developers Guide](share/doc/DEVELOPERS_GUIDE.md)
4. [VOTCA_LANGUAGE_GUIDE](share/doc/VOTCA_LANGUAGE_GUIDE.md)
5. [Code of Conduct](share/doc/CODE_OF_CONDUCT.md)

If you want to install CTP 

    prefix=WHERE/TO/INSTALL/VOTCA
    git clone -b ctp --recursive https://github.com/votca/votca.git
    cd votca
    mkdir build
    cd build
    cmake -DBUILD_CTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} ..
    cmake --build . -- -j<number of cores>
    cmake --build . --target install
