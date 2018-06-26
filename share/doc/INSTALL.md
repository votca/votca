# VOTCA Installation Guide

## VOTCA structure

VOTCA is currently composed of four major repositories:

* TOOLS - general tools 
* CSG - topographical classes and course graining functionality
* CTP - charge transport classes and functionality
* XTP - exciton and excited states functionality GW-BSE and DFT engine

TOOLS -> CSG -> CTP -> XTP

## Dependencies

Each of these repositories has different dependencies shown in the table below in alphabetical order. They are marked as following:

* r - required 
* o - optional
* n - not needed
* m - needed for building the manual
* M - needed for manpage and manual

 Dependency          | csg | ctp | xtp |
 ------------------- | --- | --- | --- |
 cmake               | r   | r   | r   |
 ghostscript         | m   | m   | m   |
 git                 | o   | o   | o   |
 graphviz            | n   | n   | r   |
 gromacs-dev         | r   | n   | n   |
 gsfonts-X11         | m   | m   | m   |
 g++                 | r   | r   | r   |
 inkscape            | n   | m   | m   |
 libboost-all-dev    | r   | r   | r   |
 libceres-dev        | n   | n   | r   | 
 libeigen3-dev       | n   | n   | r   |
 libexpat-dev        | r   | r   | r   |
 libfftw3-dev        | r   | r   | r   |
 libgsl-dev          | o   | n   | r   |
 libhdf5-dev         | o   | n   | r   |
 libsqlite3-dev      | o   | o   | o   |
 libxc-dev           | n   | n   | r   |
 pkg-config          | o   | o   | o   |
 psmisc              | r   | n   | r   |
 texlive             | m   | m   | m   |
 texlive-humanities  | m   | m   | m   |
 texlive-latex-extra | m   | m   | m   |
 txt2tags            | M   | m   | m   |
 xfig                | n   | n   | m   |
 
### Dependency Installation
#### Ubuntu
Dependencies for core functionality

    sudo apt-get install cmake git g++ libexpat-dev libfftw3-dev libgsl-dev libboost-all-dev txt2tags libsqlite3-dev libhdf5-dev graphviz pkg-config psmisc libeigen3-dev libxc-dev libceres-dev 

Dependencies for Manual

    sudo apt-get install xfig inkscape gsfonts-X11 ghostscript texlive texlive-latex-extra texlive-humanities

## Installation 

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
    
## Docker 

Votca is also available through docker and can be accessed and run with the following docker commands:

    docker pull votca/votca
    docker run -it votca/votca /bin/bash
    
    
