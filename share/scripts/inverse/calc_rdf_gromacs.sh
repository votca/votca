#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the rdf for gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo USES: get_from_mdp csg_get_interaction_property csg_get_property awk die log run_or_exit g_rdf csg_resample
   echo NEEDS: cg.inverse.gromacs.equi_time type1 type2 name step min max
   exit 0
fi

check_deps "$0"

nsteps=$(get_from_mdp nsteps) 
dt=$(get_from_mdp dt)
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time)"
equi_time=${equi_time%\%}
equi=$(awk "BEGIN{ print ${equi_time}/100*${nsteps}*${dt} }") || die "${0##*/} awk failed"

type1=$(csg_get_interaction_property type1)
type2=$(csg_get_interaction_property type2)
name=$(csg_get_interaction_property name)
binsize=$(csg_get_interaction_property step)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)

log "Running g_rdf for ${type1}-${type2}"
run_or_exit "echo -e \"${type1}\\n${type2}\" | g_rdf -b ${equi} -noxvgr -n index.ndx -bin ${binsize} -o ${name}.dist.new.xvg" -s topol.tpr
#gromacs always append xvg
run_or_exit csg_resample --in ${name}.dist.new.xvg --out ${name}.dist.new --grid ${min}:${binsize}:${max}

