#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This is a intialize stuff for gromacs
   echo Usage: ${0##*/}
   echo Needs: -
   exit 0
fi

cp ../conf.gro confout.gro