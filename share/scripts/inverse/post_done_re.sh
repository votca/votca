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
This script implements the steps required when all iterations of the relative entropy method are done.

Usage: ${0##*/}
EOF
   exit 0
fi

# get names of CG interactions
names="$(csg_get_interaction_property --all name)"

# get the last iteration directory  
step_dir="$(get_current_step_dir)"
# get the last iteration number
step_nr=$(get_step_nr $step_dir)

# copy all parameters table in the last step to the main directory
#TODO why not for_all ?
for name in $names; do
  critical cp $step_dir/${name}.param.cur .
done

# average potential parameters over these many last iteration steps
avg_steps_nr=$(csg_get_property cg.inverse.re.avg_steps '1')
is_int "${avg_steps_nr}" || die "${0##*/} cg.inverse.re.avg_steps should be an integer"
[[ ${avg_steps_nr} -gt 0 ]] || die "${0##*/} cg.inverse.re.avg_steps should be greater 0"

# sum last avg_steps_nr potential parameters
for ((i=1;i<$avg_steps_nr;i++)); do
  ((step_nr--))
  step_dir="$(get_stepname $step_nr)"
  #TODO why not for_all ?
  for name in $names; do
    do_external table combine --op + ${name}.param.cur $step_dir/${name}.param.cur ${name}.param.cur
  done
done

# now to take avg, scale final parameters table with avg_steps_nr
inv_avg_steps_nr=$(csg_calc "1.0" "/" "${avg_steps_nr}")
#TODO why not for_all ?
for name in $names; do
  do_external table linearop ${name}.param.cur ${name}.param.cur $inv_avg_steps_nr 0
done

msg --color green "Generating potential tables from the initial parameters"

#TODO csg_reupdate should get an --ext options (like csg_stat)
#     and --param-ext for parameter in-put, and --rdf-ext for rdf input
#     not many hard-coded input and output names...
critical csg_reupdate --gentable true --options $CSGXMLFILE

# rename potential tables from *.new to *.final
#TODO why not for_all ?
for name in $names; do
  critical mv ${name}.pot.new ${name}.pot.cur
  critical rm ${name}.param.new
done

