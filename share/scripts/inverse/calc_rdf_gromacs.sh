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
#20 % is warmup
equi=$(awk "BEGIN{ print 0.2*${nsteps}*${dt} }") || die "${0##*/} awk failed"

type1=$($csg_get type1)
type2=$($csg_get type2)
binsize=$($csg_get step)
log "Running g_rdf for ${type1}-${type2}"
echo -e "${type1}\n${type2}" | g_rdf -b ${equi} -n index.ndx -bin ${binsize} -o ${type1}_${type2}.dist.new.xvg >> ${CSGLOG} 2>&1 || \
die "${0##*/}: Error at g_rdf ${type1}-${type2}"
mv ${type1}_${type2}.dist.new.xvg ${type1}_${type2}.dist.new || die "${0##*/}: mv failed"

