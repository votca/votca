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
This script implemtents the pressure update

Usage: ${0##*/}

USES:  die csg_get_property do_external csg_get_interaction_property log run_or_exit cp check_deps get_current_step_nr

NEEDS: cg.inverse.program name

OPTIONAL: inverse.post_update_options.pressure.type inverse.post_update_options.pressure.do
EOF
   exit 0
fi

check_deps "$0"

step_nr="$(get_current_step_nr)"
sim_prog="$(csg_get_property cg.inverse.program)"
name=$(csg_get_interaction_property name)

p_now="$(do_external pressure $sim_prog)" || die "${0##*/}: do_external pressure $sim_prog failed"
[ -z "$p_now" ] && die "${0##*/}: Could not get pressure from simulation"
log "New pressure $p_now"

ptype="$(csg_get_interaction_property inverse.post_update_options.pressure.type simple)"
pscheme=( $(csg_get_interaction_property inverse.post_update_options.pressure.do 1 ) )
pscheme_nr=$(( ( $step_nr - 1 ) % ${#pscheme[@]} ))

if [ "${pscheme[$pscheme_nr]}" = 1 ]; then
   log "Apply ${ptype} pressure correction for interaction ${name}"
   do_external pressure_cor $ptype $p_now pressure_cor.d
   do_external table add pressure_cor.d ${name}.dpot.cur ${name}.dpot.new
else
   log "NO pressure correction for interaction ${name}"
   run_or_exit cp ${name}.dpot.cur ${name}.dpot.new
fi
