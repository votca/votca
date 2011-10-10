#! /bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script calcs the thermoforce out of gromacs density for the AdResS therm force

Usage: ${0##*/} infile outfile
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

[[ -f $2 ]] && die "${0##*/}: $2 is already there"

name=$(csg_get_interaction_property name)

mdp="$(csg_get_property cg.inverse.gromacs.mdp)"
#this needs to done in step_000
[ -f "$mdp" ] || cp_from_main_dir "$mdp"

infile="${1}"
endfile="${2}"

adress_type=$(get_simulation_setting adress_type)
if [ $adress_type = "sphere" ]; then
  #note: in the spehere case (no symmetrizing necessary) infile stays dens.${name}.xvg, so this gets used for next step
  :
else
  outfile="${name}.sym.dens"
  adressc="$(get_simulation_setting adress_reference_coords "0")"
  ref="$(echo "$adressc" | awk '{if (NF<1) exit 1; print $1;}')" || die "${0##*/}: we need at least one number in adress_reference_coords, but got '$adressc'"
  critical do_external density symmetrize --infile "$infile" --outfile "$outfile" --adressc "$ref"
  infile="${outfile}"
fi

comment="$(get_table_comment)"

min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step="$(csg_get_interaction_property step)"

sp_min=$(csg_get_interaction_property inverse.tf.spline_start )
sp_max=$(csg_get_interaction_property inverse.tf.spline_end )
sp_step="$(csg_get_interaction_property inverse.tf.spline_step)"

#resample to a bigger grid
bigger="${name}.extended.dens"
critical csg_resample --type cubic --in "$infile" --out "$bigger" --grid "$sp_min:$step:$sp_max" --comment "$comment"

#calculate derivative of the density using csg_resample on a spline grid
forcefile="${name}.tf"
smooth="${name}.smooth.dens"
critical csg_resample --type cubic --in "$bigger" --out "$smooth" --grid "$sp_min:$step:$sp_max" --derivative "$forcefile" --fitgrid "$sp_min:$sp_step:$sp_max" --comment "$comment"

#multiply the prefactor on
prefactor="$(csg_get_interaction_property inverse.tf.prefactor)"
forcefile_pref="${name}.tf_with_prefactor"
cg_prefactor="$(csg_get_interaction_property --allow-empty inverse.tf.cg_prefactor)"
if [[ -z $cg_prefactor ]]; then 
  echo "Using fixed prefactor $prefactor" 
  do_external table linearop "$forcefile" "${forcefile_pref}" "${prefactor}" 0.0
else
  echo "Using linear interpolation of prefactors between $prefactor and $cg_prefactor"
  do_external table scale "$forcefile" "$forcefile_pref" "$prefactor" "$cg_prefactor"
fi

#cut it down to the range min to max
forcefile_smooth="${name}.tf_smooth"
do_external table smooth_borders --infile "$forcefile_pref" --outfile "$forcefile_smooth" --xstart "$min" --xstop "$max"

#integrate the force table
do_external table integrate "$forcefile_smooth" "${endfile}"
#we do not need to multiply it with -1, because we assue that the prefactor is always positiv, but the update is -prefacor * grad density

