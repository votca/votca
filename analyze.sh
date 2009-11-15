#!/bin/bash

msg "updating plots"

name=$(csg_get_interaction_property name)
cp $name.pot.cur $name.pot.new

cp $name.dist.new ..
cp $name.dpot.new ..
cp $name.pot.new ..

dir=$PWD
cd ..
./convergence.sh > convergence.dat 2>/dev/null
killall -9 gnuplot_x11 &> /dev/null
./plot_all.sh
cd $dir

