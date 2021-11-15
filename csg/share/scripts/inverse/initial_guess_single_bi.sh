#! /bin/bash
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

show_help () {
  cat <<EOF
${0##*/}, version %version%
This script prepares a single interaction with Boltzmann inversion

Usage: ${0##*/} [--help]
EOF
}

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step=$(csg_get_interaction_property step )
comment="$(get_table_comment)"
bondtype="$(csg_get_interaction_property bondtype)"
output="${name}.pot.new"

target=$(csg_get_interaction_property inverse.target)
msg "Using initial guess from dist ${target} for ${name}"
# resample target dist
do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt" 
# initial guess from rdf
raw="$(critical mktemp ${name}.pot.new.raw.XXX)"
multistate="$(csg_get_property cg.inverse.multistate.enabled)"
if [[ $multistate == true ]]; then
  state_nr=$(get_state_nr)
  kbt_array=( $(csg_get_property cg.inverse.multistate.state_kBTs) )
  kbt="${kbt_array[state_nr]}"
else
  kbt="$(csg_get_property cg.inverse.kBT)"
fi
dist_min="$(csg_get_property cg.inverse.dist_min)"
do_external dist invert --type "${bondtype}" --kbT "${kbt}" --min "${dist_min}" ${name}.dist.tgt ${raw}
# smooth
smooth="$(critical mktemp ${name}.pot.new.smooth.XXX)"
critical csg_resample --in ${raw} --out ${smooth} --grid ${min}:${step}:${max} --comment "${comment}"
# extrapolate
extrapolate="$(critical mktemp ${name}.pot.new.extrapolate.XXX)"
do_external potential extrapolate --type "$bondtype" "${smooth}" "${extrapolate}"
# shift
shifted="$(critical mktemp ${name}.pot.new.shifted.XXX)"
do_external potential shift --type "${bondtype}" ${extrapolate} ${shifted}
# scale new potentials
if [[ ${bondtype} == "non-bonded" ]]; then
  scale_param="scale_non_bonded"
else
  scale_param="scale_bonded"
fi
scaling_factor="$(csg_get_property "cg.inverse.initial_guess.${scale_param}")"
if $(awk "BEGIN {exit (${scaling_factor} != 1.0 ? 0 : 1)}"); then
  scaled="$(critical mktemp ${name}.pot.new.scaled.XXX)"
  do_external table linearop "${shifted}" "${scaled}" "${scaling_factor}" "0"
else
  scaled="${shifted}"
fi
# set all flags to 'i'
do_external table change_flag "${scaled}" "${output}"
