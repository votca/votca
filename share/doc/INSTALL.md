# VOTCA Installation Guide

## VOTCA structure

VOTCA is currently composed of four major repositories:

* TOOLS - general tools 
* CSG - topographical classes and course graining functionality
* CTP - charge transport classes and functionality
* XTP - exciton and excited states functionality GW-BSE and DFT engine

TOOLS -> CSG -> CTP -> XTP

### Dependency Installation
#### Ubuntu
Dependencies for core functionality

    sudo apt-get install cmake git g++ libexpat-dev libfftw3-dev libgsl-dev libboost-all-dev txt2tags libsqlite3-dev libhdf5-dev graphviz pkg-config psmisc libeigen3-dev libxc-dev libceres-dev 

Dependencies for Manual

    sudo apt-get install xfig inkscape gsfonts-X11 ghostscript texlive texlive-latex-extra texlive-humanities

## General Installation Instructions 

To install the full package:

    prefix=WHERE/TO/INSTALL/VOTCA
    git clone --recursive https://github.com/votca/votca.git
    cd votca
    mkdir build
    cd build
    cmake -DBUILD_CSGAPPS=ON -DBUILD_CTP=ON -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} ..
    make -j5

### cmake Flags

* `BUILD_CSGAPPS` - Build the extra csg applications repo (ON/OFF, Default OFF)
* `BUILD_XTP` - Build the xtp repo (ON/OFF, Default OFF)
* `BUILD_CTP` - Build the ctp repo (ON/OFF, Default OFF)
* `CMAKE_INSTALL_PREFIX` - where to install the votca executables (Default is /usr/local/bin)
* `ENABLE_TESTING` - compile tests (ON/OFF, Default OFF)

## Yum (Native - CentOS, Fedora)

   yum install epel-release
   yum update
   yum install votca-csg
   yum install votca-xtp

## Zypper (Native - OpenSuse, SLES)

   zypper install votca-csg
   zypper install votca-xtp

## Apt-get (Native - Debian, Ubuntu)

    apt-get install votca-csg
    
## Emerge (Native - Gentoo)

    emerge votca-csg

## Spack

    git clone clone https://github.com/spack/spack.git spack
    source spack/share/spack/setup-env.sh
    spack install votca-csg
    spack install votca-xtp

## Docker 

Votca is also available through docker and can be accessed and run with the following docker commands:

    docker pull votca/votca
    docker run -it votca/votca /bin/bash
    
    
