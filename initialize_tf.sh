#! /bin/bash

if [ "$1" = "--help" ]; then
  echo This script implemtents the function initialize
  echo for the Inverse Boltzmann Method
  echo Usage: ${0##*/}
  echo USES: do_external csg_get_interaction_property log run_or_exit csg_resample log
  echo NEEDS: name min max step
  exit 0
fi

check_deps "$0"

 
name=$(csg_get_interaction_property name)
#dummy pot

 if [ -f ../${name}.pot.in ]; then
   log "Using given thermforce ${name}.pot.in for ${name}"
   min=$(csg_get_interaction_property min )
   max=$(csg_get_interaction_property max )
   step=$(csg_get_interaction_property step )
   run_or_exit csg_resample --in ../${name}.pot.in --out ${name}.pot.new --grid ${min}:${step}:${max}
 elif [ -f ../dens_${name}.xvg ]; then
   log "Calculating initial therm force from input density file dens_${name}.xvg"
    xstart="$(csg_get_property cg.tf.xstart)"
    xstop="$(csg_get_property cg.tf.xstop)"
    step="$(csg_get_interaction_property step)"
    prefactor="$(csg_get_property cg.tf.prefactor)"
    splinestep="$(csg_get_property cg.tf.splinesmoothstep)"
    
    adressw="$(csg_get_property cg.tf.adressw)"
    adressh="$(csg_get_property cg.tf.adressh)"
    adressc="$(csg_get_property cg.tf.adressc)"
    
    infile="dens_$name.xvg"
    outfile="dens_$name\_smooth.xvg"
    forcefile="thermforce_$name.xvg"
    forcefile_smooth="thermforce_smooth_$name.xvg"

    run_or_exit csg_resample --in ../$infile --out $outfile --grid $xstart:$step:$xstop --derivative $forcefile --spfit $xstart:$splinestep:$xstop
    
    run_or_exit do_external table smooth_borders --infile $forcefile --outfile $forcefile_smooth --xstart $xstart --xstop $xstop
    
    run_or_exit do_external table integrate $forcefile_smooth ${name}.pot.new $prefactor
 else
	die "initial therm_force_file ${name}.pot.in not found or initial density file dens_${name}.xvg not found" 
   # RDF_to_POT.pl just does log g(r) + extrapolation
   # log "Using intial guess from RDF for ${name}"
   #run_or_exit do_external rdf pot ${name}.dist.tgt ${name}.pot.new
 fi

