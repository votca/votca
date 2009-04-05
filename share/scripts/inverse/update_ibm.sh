#! /bin/bash

if [ "$1" = "--help" ]; then
   echo This script implemtents the function update
   echo for the Inverse Boltzmann Method
   echo Usage: ${0##*/} step_nr
   echo Needs:  run_or_exit, \$source_wrapper, update_POT.pl
   exit 0
fi

if [ -z "$1" ]; then
  echo Missing argument for ${0##*/} > /dev/stderr
  exit 1
fi

if [ -z "$SOURCE_WRAPPER" ]; then
   echo source_wrapper not defined > /dev/stderr
   exit 1
fi

if [ $(type -t run_or_exit) != "function" ]; then
   echo Could not find function run_or_exit > /dev/stderr
   exit 1
fi

update_POT="$($SOURCE_WRAPPER --direct update_POT.pl)" || exit 1
add_POT=$($SOURCE_WRAPPER --direct add_POT.pl) || exit 1

csg="csg_property --file $CSGXMLFILE --short --path cg.${bondtype} --filter name=${name}"
scheme=( $($csg --print .ibm.do_potential | sed '1d' ) )
pscheme=( $($csg --print .ibm.do_pressure | sed '1d' ) )
scheme_nr=$(( ( $1 - 1 ) % ${#scheme[@]} ))
pscheme_nr=$(( ( $1 - 1 ) % ${#pscheme[@]} ))

if [ "${scheme[$scheme_nr]}" = 1 ]; then
   echo Update potential ${type1}-${type2} : yes 
   #update ibm
   run_or_exit --log log_update_POT_${type1}_${type2} $update_POT rdf_${type1}_${type2}_aim.xvg rdf_${type1}_${type2}.xvg delta_pot_${type1}_${type2}.d
   run_or_exit --log log_add_POT_${type1}_${type2} $add_POT table_${type1}_${type2}.d delta_pot_${type1}_${type2}.d table_${type1}_${type2}_new1.d
else
   echo Update potential ${type1}-${type2} : no 
   cp table_${type1}_${type2}.d table_${type1}_${type2}_new1.d
fi
if [ "${pscheme[$pscheme_nr]}" = 1 ]; then
   echo Update presuure ${type1}-${type2} : yes
   run_or_exit --log log_add_POT2_${type1}_${type2} $add_POT table_${type1}_${type2}_new1.d pressure_cor.d table_${type1}_${type2}_new.d
else
   echo Update pressure ${type1}-${type2} : no 
   cp  table_${type1}_${type2}_new1.d table_${type1}_${type2}_new.d
fi
