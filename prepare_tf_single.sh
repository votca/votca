#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This script implemtents the function initialize
  echo for the Inverse Boltzmann Method
  echo Usage: ${0##*/}
  echo USES: do_external csg_get_interaction_property run_or_exit csg_resample
  echo NEEDS: name min max step
  exit 0
fi

check_deps "$0"


r_cut=$(csg_get_interaction_property max)
adressh="$(csg_get_property cg.tf.adressh)"

[ $(echo $r_cut == $adressh|bc) -eq 0 ] && die "r_cut $r_cut is not equal adressh $adressh! This has to be set manually"

 
name=$(csg_get_interaction_property name)
#dummy pot

 if [ -f ../${name}.pot.in ]; then
   echo "Using given thermforce ${name}.pot.in for ${name}"
   min=$(csg_get_interaction_property min )
   max=$(csg_get_interaction_property max )
   step=$(csg_get_interaction_property step )
   comment="$(get_table_comment)"
   run_or_exit csg_resample --in ../${name}.pot.in --out ${name}.pot.new --grid ${min}:${step}:${max} --comment "$comment"
 elif [ -f ../dens.${name}.xvg ]; then
   echo "Calculating initial therm force from input density file dens_${name}.xvg"
   cp ../dens.$name.xvg .
   run_or_exit do_external calc thermforce
   cp ${name}.dpot.new ${name}.pot.new 
 else
	die "initial therm_force_file ${name}.pot.in not found or initial density file dens.${name}.xvg not found" 
   # RDF_to_POT.pl just does log  g(r) + extrapolation
   # log "Using intial guess from RDF for ${name}"
   #run_or_exit do_external rdf pot ${name}.dist.tgt ${name}.pot.new
 fi

