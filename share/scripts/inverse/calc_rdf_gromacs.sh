#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the rdf for gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo Needs: g_rdf
   exit 0
fi

for exe in g_rdf awk; do
   [[ -n $(type -p $exe) ]] || die "${0##*/}: Could not find $exe"
done

nsteps=$(get_from_mdp nsteps) 
dt=$(get_from_mdp dt)
equi_time="$(get_sim_property gromacs.equi_time)"
equi_time=${equi_time%\%}
equi=$(awk "BEGIN{ print ${equi_time}/100*${nsteps}*${dt} }") || die "${0##*/} awk failed"

type1=$($csg_get type1)
type2=$($csg_get type2)
binsize=$($csg_get step)
log "Running g_rdf for ${type1}-${type2}"
run_or_exit "echo -e \"${type1}\\n${type2}\" | g_rdf -b ${equi} -n index.ndx -bin ${binsize} -o ${type1}_${type2}.dist.new.xvg" 
#gromacs allways append xvg
mv ${type1}_${type2}.dist.new.xvg ${type1}_${type2}.dist.new || die "${0##*/}: mv failed"

