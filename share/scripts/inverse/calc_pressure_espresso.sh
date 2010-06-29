#! /bin/bash
#
# Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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
This script calcs the pressure for espresso
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES: log check_deps

NEEDS: 
EOF
   exit 0
fi

#####################
# NOT IMPLEMENTED YET
#####################
log "!! Pressure correction in ESPResSo not yet implemented !!"

check_deps "$0"



# mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
# [ -f "$mdp" ] || die "${0##*/}: gromacs mdp file '$mdp' not found"

# tpr="$(csg_get_property cg.inverse.gromacs.g_energy.topol "topol.tpr")"
# [ -f "$tpr" ] || die "${0##*/}: Gromacs tpr file '$tpr' not found"

# opts="$(csg_get_property --allow-empty cg.inverse.gromacs.g_energy.opts)"

# nsteps=$(get_from_mdp nsteps "$mdp")
# dt=$(get_from_mdp dt "$mdp")
# equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
# first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"

# begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"

# log "Running g_energy"
# echo Pressure | run_or_exit g_energy -b "${begin}" -s "${tpr}" ${opts}
# #the number pattern '-\?[0-9][^[:space:]]*[0-9]' is ugly, but it supports X X.X X.Xe+X Xe-X and so on
# p_now=$(csg_taillog -30 | sed -n 's/^Pressure[^-0-9]*\(-\?[0-9][^[:space:]]*[0-9]\)[[:space:]].*$/\1/p' ) || \
#   die "${0##*/}: awk failed"
# [ -z "$p_now" ] && die "${0##*/}: Could not get pressure from simulation"
# echo ${p_now}
