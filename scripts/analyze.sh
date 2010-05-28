#!/bin/bash


cp $1 $2

name=$(csg_get_interaction_property name)
cp $name.dpot.new ..

cd $(get_main_dir)
do_external postadd convergence > convergence.dat 2>/dev/null
msg "updating plots"
killall -9 gnuplot_x11 &> /dev/null
do_external postadd plot

