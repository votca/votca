[![Codacy Badge](https://api.codacy.com/project/badge/Grade/9f3db8fd6bd243128b42788a308250ae)](https://app.codacy.com/app/JoshuaSBrown/votca?utm_source=github.com&utm_medium=referral&utm_content=votca/votca&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.org/votca/votca.svg?branch=master)](https://travis-ci.org/votca/votca)
[![pipeline status](https://gitlab.com/votca/votca/badges/master/pipeline.svg)](https://gitlab.com/votca/votca/commits/master)

This is VOTCA's next generation build system.

Usage:

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'stable', 'ctp' or 'v1.4.1'
    git clone -b ${version} --recursive https://github.com/votca/votca.git
    cd votca
    mkdir build
    cd build
    cmake -DBUILD_CSGAPPS=ON -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} ..
    make -j<number of cores>
    make install

Using this code via docker:

    docker pull votca/votca
    docker run -it votca/votca /bin/bash

For further details see:

1. [Installation](share/doc/INSTALL.md)
2. [Further Information](http://www.votca.org)
3. [Developers Guide](share/doc/DEVELOPERS_GUIDE.md)
4. [Code of Conduct](share/doc/CODE_OF_CONDUCT.md)
