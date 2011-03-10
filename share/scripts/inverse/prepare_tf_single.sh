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
This script implemtents the function initialize for the thermodynamic force iteration

Usage: ${0##*/}
EOF
exit
fi

r_cut="$(csg_get_interaction_property max)"
#TODO compare this to mdp file....
adressh="$(csg_get_property cg.tf.adressh)"

res="$(awk -v x="$r_cut" -v y="$adressh" 'BEGIN{ print (sqrt((x-y)**2)>0.001)?1:0}')" || die "${0##*/}: awk failed"
[ "$res" != "0" ] && die "Error: r_cut $r_cut is not equal adressh $adressh"

name=$(csg_get_interaction_property name)

if [ -f ../${name}.pot.in ]; then
   msg "Using given thermforce ${name}.pot.in for ${name}"
   min=$(csg_get_interaction_property min )
   max=$(csg_get_interaction_property max )
   step=$(csg_get_interaction_property step )
   comment="$(get_table_comment)"
   critical csg_resample --in ../${name}.pot.in --out ${name}.pot.new --grid ${min}:${step}:${max} --comment "$comment"
elif [ -f ../dens.${name}.xvg ]; then
   msg "Calculating initial therm force from input density file dens_${name}.xvg"
   critical cp ../dens.$name.xvg .
   critical do_external calc thermforce
   #TODO call csg_resample here
   critical cp ${name}.dpot.new ${name}.pot.new 
else
   die "initial therm_force_file ${name}.pot.in or initial density file dens.${name}.xvg not found" 
fi

