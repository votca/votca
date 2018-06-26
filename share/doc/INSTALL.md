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

    apt-get install make cmake git g++ libexpat-dev libfftw3-dev libgsl-dev libboost-all-dev txt2tags libsqlite3-dev octave gnuplot python-numpy libhdf5-dev graphviz pkg-config psmisc libint2-dev libeigen3-dev libxc-dev libceres-dev
    
Dependencies for Manual

    apt-get install ghostscript texlive inkscape transfig texlive-latex-extra texlive-pstricks
    
### Fedora
Dependencies for core functionality

     dnf install make cmake git gcc-c++ expat-devel fftw-devel gsl-devel boost-devel txt2tags sqlite-devel procps-ng octave gnuplot python2-numpy psmisc hdf5-devel lammps libint2-devel eigen3-devel libxc-devel ceres-solver-devel python-numpy

Dependencies for Manual

     dnf install ghostscript texlive doxygen texlive-appendix texlive-wrapfig texlive-a4wide texlive-xstring inkscape transfig texlive-units texlive-sidecap texlive-bclogo texlive-mdframed texlive-braket graphviz ImageMagick ghostscript-tools-dvipdf


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

## Fedora

    dnf install votca-csg votca-xtp

## CentOs

    yum install epel-release
    yum update
    yum install votca-csg

## OpenSuse

    zypper install votca-csg votca-xtp
    
## SLES

    SUSEConnect -p PackageHub/12.2/x86_64
    zypper install votca-csg    

## Debian / Ubuntu

    apt-get install votca-csg
    
## Gentoo 

    emerge votca-csg votca-xtp

## Spack

    git clone clone https://github.com/spack/spack.git spack
    source spack/share/spack/setup-env.sh
    spack install votca-csg
    spack install votca-xtp

## Docker 

Votca is also available through docker and can be accessed and run with the following docker commands:

    docker pull votca/votca
    docker run -it votca/votca /bin/bash
    
    
