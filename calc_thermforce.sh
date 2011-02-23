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

Usage: ${0##*/}
EOF
   exit 0
fi

name=$(csg_get_interaction_property name)

splinedelta="$(csg_get_property cg.tf.splinedelta)"
step="$(csg_get_interaction_property step)"
prefactor="$(csg_get_property cg.tf.prefactor)"

cg_prefactor="$(csg_get_property --allow-empty cg.tf.cg_prefactor)"
splinestep="$(csg_get_property cg.tf.splinesmoothstep)"

#TODO compare this to mdp file
adressw="$(csg_get_property cg.tf.adressw)"
adressh="$(csg_get_property cg.tf.adressh)"
adressc="$(csg_get_property cg.tf.adressc)"
adress_type=$(get_from_mdp adress_type "$mdp")

mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
#this needs to done in step_000
if [ ! -f "$mdp" ]; then
  cp_from_maindir "$mdp"
fi

infile="dens.${name}.xvg"
outfile="dens.${name}.smooth.xvg"

xstart="$(awk -v x="$adressc" -v y="$adressw" 'BEGIN{print x+y}')" || die "${0##*/}: awk failed"
xstop="$(awk -v x="$xstart" -v y="$adressh" 'BEGIN{print x+y}')" || die "${0##*/}: awk failed"

infile="dens.${name}.xvg"
if [ $adress_type = "sphere" ]; then
  #note : in the spehere case (no symmetrizing necessary) infile stays dens.${name}.xvg, so this gets used for next step
  :
else
  outfile="dens.${name}.symm.xvg"
  critical do_external density symmetrize --infile "$infile" --outfile "$outfile" --adressc "$adressc"
  infile="dens.${name}.symm.xvg"
fi
outfile="dens.${name}.smooth.xvg"

#TODO remove this hardcoded stuff
forcefile="thermforce.${name}.xvg"
forcefile_pref="thermforce.wpref.${name}.xvg"
forcefile_smooth="thermforce.smooth.${name}.xvg"

spxstart="$(awk -v x="$xstart" -v y="$splinedelta" 'BEGIN{print x-y}')" || die "${0##*/}: awk failed"
spxstop="$(awk -v x="$xstop" -v y="$splinedelta" 'BEGIN{print x+y}')" || die "${0##*/}: awk failed"

comment="$(get_table_comment)"

critical csg_resample --type cubic --in "$infile" --out "$outfile" --grid "$spxstart:$step:$spxstop" --derivative "$forcefile" --fitgrid "$spxstart:$splinestep:$spxstop" --comment "$comment"

if [ -z "$cg_prefactor" ];then
       echo "Using fixed prefactor $prefactor "	
       critical do_external tf apply_prefactor $forcefile $forcefile_pref $prefactor
else
       echo "Using linear interpolation of prefactors. Ex. pref: $prefactor CG. pref : $cg_prefactor"
       critical do_external tf apply_prefactor $forcefile $forcefile_pref $prefactor $cg_prefactor
fi

critical do_external table smooth_borders --infile "$forcefile_pref" --outfile "$forcefile_smooth" --xstart "$xstart" --xstop "$xstop"

critical do_external table integrate "$forcefile_smooth" "minus_${name}.dpot.new"
critical do_external table linearop "minus_${name}.dpot.new" "minus_${name}.dpot.new" -1.0 0.0

