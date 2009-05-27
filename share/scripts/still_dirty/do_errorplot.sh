#!/bin/bash

do_errors() {
  csg_imc --top topol.tpr --trj traj.xtc --cg ../water_cg.xml --options ../cg.xml \
    --do-imc --begin 100 \
    --write-every $1 --do-blocks --nframes $((16*$1))
  ~/src/csg/share/scripts/inverse/errorbars.sh ibm
  cp CG-CG.dpot.err ibm.err.$i
  ~/src/csg/share/scripts/inverse/errorbars.sh imc
  cp CG-CG.dpot.err imc.err.$i
}

for i in $(seq 1 250); do
  do_errors $i
done
