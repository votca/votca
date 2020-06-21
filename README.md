[![Codacy Badge](https://api.codacy.com/project/badge/Grade/48a26be8dd8b4f0fa67c93646fa6d30d)](https://www.codacy.com/manual/votca-package/votca?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=votca/votca&amp;utm_campaign=Badge_Grade)
[![CI](https://github.com/votca/votca/workflows/CI/badge.svg)](https://github.com/votca/votca/actions?query=branch%3Agithub_actions+workflow%3ACI)
[![DOI](https://zenodo.org/badge/75022030.svg)](https://zenodo.org/badge/latestdoi/75022030)

This is VOTCA's next generation build system for CSG and XTP. It allows you to easily install: 

-   VOTCA-CSG, a library which provides tools to develop coarse-grained potentials from atomistic simulation data
-   VOTCA-XTP, a library designed to determine electronic properties of organic materials from atomistic MD-trajectories.

Usage:

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'stable' or 'v1.6.1'
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

You find the VOTCA-CTP repository [here](https://gitlab.mpcdf.mpg.de/votca/votca)
