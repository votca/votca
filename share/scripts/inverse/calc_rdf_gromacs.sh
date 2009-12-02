#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script calcs the rdf for gromacs
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/}
   echo USES: get_from_mdp csg_get_interaction_property csg_get_property awk log run_or_exit g_rdf csg_resample is_done mark_done msg
   echo NEEDS: type1 type2 name step min max
   echo OPTIONAL: cg.inverse.gromacs.equi_time cg.inverse.gromacs.first_frame 
   exit 0
fi

check_deps "$0"

dt=$(get_from_mdp dt)
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"

type1=$(csg_get_interaction_property type1)
type2=$(csg_get_interaction_property type2)
name=$(csg_get_interaction_property name)
binsize=$(csg_get_interaction_property step)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)

begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"

log "Running g_rdf for ${type1}-${type2}"
if is_done "rdf-$name"; then
  msg "g_rdf for ${type1}-${type2} is already done"
else
  run_or_exit "echo -e \"${type1}\\n${type2}\" | g_rdf -b ${begin} -noxvgr -n index.ndx -bin ${binsize} -o ${name}.dist.new.xvg -s topol.tpr"
#gromacs always append xvg
  run_or_exit csg_resample --in ${name}.dist.new.xvg --out ${name}.dist.new --grid ${min}:${binsize}:${max}
  mark_done "rdf-$name"
fi

