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

    apt-get install make cmake git g++ libexpat-dev libfftw3-dev libgsl-dev libboost-all-dev txt2tags libsqlite3-dev octave gnuplot python-numpy libhdf5-dev graphviz pkg-config psmisc libint2-dev libeigen3-dev libxc-dev libceres-dev libgromacs-dev gromacs-openmpi
    
Dependencies for Manual

    apt-get install ghostscript texlive inkscape transfig texlive-latex-extra texlive-pstricks
    
#### Fedora
Dependencies for core functionality

     dnf install make cmake git gcc-c++ expat-devel fftw-devel gsl-devel boost-devel txt2tags sqlite-devel procps-ng octave gnuplot python2-numpy psmisc hdf5-devel lammps libint2-devel eigen3-devel libxc-devel ceres-solver-devel python-numpy gromacs-devel gromacs gromacs-openmpi

Dependencies for Manual

     dnf install ghostscript texlive doxygen texlive-appendix texlive-wrapfig texlive-a4wide texlive-xstring inkscape transfig texlive-units texlive-sidecap texlive-bclogo texlive-mdframed texlive-braket graphviz ImageMagick ghostscript-tools-dvipdf


## General (Source) Installation Instructions 

To install the full package:

    prefix=WHERE/TO/INSTALL/VOTCA
    version=master # or 'stable' or 'v1.4.1'
    git clone -b ${version} --recursive https://github.com/votca/votca.git
    cd votca
    mkdir build
    cd build
    cmake -DBUILD_CSGAPPS=ON -DBUILD_CTP=ON -DBUILD_XTP=ON -DCMAKE_INSTALL_PREFIX=${prefix} ..
    make -j5
    make install

### Common CMake Flags

* `BUILD_CSGAPPS` - Build the extra csg applications repo (ON/OFF, Default OFF)
* `BUILD_XTP` - Build the xtp repo (ON/OFF, Default OFF)
* `BUILD_CTP` - Build the ctp repo (ON/OFF, Default OFF)
* `CMAKE_INSTALL_PREFIX` - where to install the votca executables (Default is /usr/local/bin)
* `ENABLE_TESTING` - compile tests (ON/OFF, Default OFF)

### Other CMake Flags

 * `BUILD_CSG_MANUAL` - Build csg pdf manual
 * `BUILD_CTP_MANUAL` - Build ctp pdf manual
 * `BUILD_XTP_MANUAL` - Build xtp pdf manual
 * `WITH_GMX` - Build with Gromacs support (ON/OFF, Default ON)

## Legacy (Source) Installation Instructions

Check [dependencies](#dependency-installation) first. Do **NOT** download anything yourself, this is done by `build.sh` below.

When installing for the first time, run the following commands in bash:
```
prefix=WHERE/TO/INSTALL/VOTCA
mkdir -p ${prefix}/src
cd ${prefix}/src
wget https://raw.githubusercontent.com/votca/buildutil/master/build.sh
chmod +x build.sh
./build.sh --prefix ${prefix} --dev tools csg ctp xtp
```
Replace the _WHERE/TO/INSTALL/VOTCA_ with directory where you want to install votca, usually ${HOME}/votca is good choice.

This commands will:
  * Create a directory, where the sources are to be stored (above `${prefix}/src`)
  * Go there
  * Download our build script
  * Make the build script executable
  * Run the build script to download all votca modules using git, it uses the development version, as xtp still undergoes extensive development and no releases are available yet.

## Binary Packages for various Linux Distributions
### Fedora

    dnf install votca-csg votca-xtp

### CentOs

    yum install epel-release
    yum update
    yum install votca-csg

### OpenSuse

    zypper install votca-csg votca-xtp
    
### SLES

    SUSEConnect -p PackageHub/12.2/x86_64
    zypper install votca-csg    

### Debian / Ubuntu

    apt-get install votca-csg
    
### Gentoo 

    emerge votca-csg votca-xtp

### Spack
[Spack](https://spack.io/) is an package manager, which allows to build VOTCA and all its dependencies:

    git clone clone https://github.com/spack/spack.git spack
    source spack/share/spack/setup-env.sh
    spack install votca-csg
    spack install votca-xtp

#### Development version
Spack can also install the latest development version from git using:

    spack install votca-csg@develop

### Docker 

Votca is also available through docker and can be accessed and run with the following docker commands:

    docker pull votca/votca
    docker run -it votca/votca /bin/bash
    
    
