#! /usr/bin/env bash
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
postadd average script,
calcs averages of (\${name}.DIST.cur) for the past few steps
and saves it to \${name}.DIST.avg
DIST can be specified by average.what option

usage: ${0##*/}
EOF
   exit 0
fi

name=$(csg_get_interaction_property name)

# these are arrays
what_to_do_list=( $(csg_get_interaction_property inverse.post_add_options.average.what) )
method="$(csg_get_property cg.inverse.method)"

# get current step directory name and number
step_dir="$(get_current_step_dir)"
step_nr=$(get_step_nr $step_dir)

# maximum number of steps to be averaged for
max_steps_nr=$(csg_get_property cg.inverse.average.steps '2')

# we allow multiple things per interaction to be averaged
for ((i=0;i<${#what_to_do_list[@]};i++)); do
  dist=${what_to_do_list[$i]}
  # pot in case of re is special
  # do not compute potential averages directly
  # compute it from the avg parameters
  if [[ ${method} = "re" && ${dist} = "pot" ]]; then
    # for dist = pot check if avg parameters have been computed or not
    if [[ -f ${name}.param.avg ]]; then
      critical csg_reupdate --gentable true --interaction "${name}" --param-in-ext param.avg --pot-out-ext pot.avg --options $CSGXMLFILE
    else
      die "${0##*/}: file '${name}.param.avg' was not found. Make sure 'param' is specified before 'pot' in the what-do list at '$name.inverse.post_add_options.average.what'."
    fi
  else
    # store the tables from last max_steps_nr steps
    for ((step_i=0;step_i<$max_steps_nr && $step_nr>=0;step_i++)); do
      step_dir="$(get_stepname $step_nr)"
      if [[ -f $(get_main_dir)/$step_dir/${name}.${dist}.cur ]]; then
        tables[$step_i]="$(get_main_dir)/$step_dir/${name}.${dist}.cur"
      fi
      ((step_nr--))
    done
    # compute the average if more than one tables found
    if [[ ${#tables[@]} -gt 1 ]]; then
      do_external table average --output ${name}.${dist}.avg "${tables[@]}"
    else
      # copy the single table to *.avg
      critical cp ${tables[0]} ${name}.${dist}.avg
    fi
  fi
done
