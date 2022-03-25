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

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script does a Newton update of non-bonded interactions.

It expects to find jacobian.npz already generated

Usage: ${0##*/}  [--help]
EOF
   exit 0
fi

sim_prog="$(csg_get_property cg.inverse.program)"

# topology for molecular conections and volume
topol=$(csg_get_property cg.inverse.topol_xml)
[[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"

# volume
volume=$(critical csg_dump --top "$topol" | grep 'Volume' | awk '{print $2}')
([[ -n "$volume" ]] && is_num "$volume") || die "could not determine the volume from file ${topol}"

# verbose
verbose=$(csg_get_property cg.inverse.verbose)
step_nr=$(get_current_step_nr)
[[ "${verbose}" == 'true' ]] && verbose_flag="--verbose"
[[ "${verbose}" == 'step0+1' ]] && [[ $step_nr == '0' || $step_nr == '1' ]] && verbose_flag="--verbose"

# Some arguments (cut_off, kBT) will be read directly from the settings.xml. They do not have a default in csg_defaults.xml.
# Others (closure, ...) could also be read from the settings file, but this bash script handles the defaults.
do_external update newton_py newton \
  ${verbose_flag-} \
  --jacobian "jacobian.npz" \
  --volume "$volume" \
  --topol "$topol" \
  --options "$CSGXMLFILE" \
  --g-tgt-ext "dist.tgt" \
  --g-cur-ext "dist.new" \
  --out "dpot.pure_n"

# resample potentials. This is needed because non-bonded.max is sometimes larger than iie.cut-off and the former should define the end of the table
for_all "non-bonded" 'csg_resample --in $(csg_get_interaction_property name).dpot.pure_n --out $(csg_get_interaction_property name).dpot.grid_adapted --grid $(csg_get_interaction_property min):$(csg_get_interaction_property step):$(csg_get_interaction_property max) --comment "adapted to grid in update_iie.sh"'
# csg_resample alone will sometimes lead to non-zero values in the tail, table extrapolate will make it zero
for_all "non-bonded" 'do_external table extrapolate --function constant --region right --no-flagupdate $(csg_get_interaction_property name).dpot.grid_adapted $(csg_get_interaction_property name).dpot.new'

# overwrite with zeros if do_potential=0
do_potential_zero_overwrite() {
  step_nr=$(get_current_step_nr)
  scheme=( $(csg_get_interaction_property inverse.do_potential) )
  scheme_nr=$(( ( $step_nr - 1 ) % ${#scheme[@]} ))
  name=$(csg_get_interaction_property name)
  if [[ ${scheme[$scheme_nr]} == 0 ]]; then
    echo "Update potential ${name} : no"
    min=$(csg_get_interaction_property min)
    max=$(csg_get_interaction_property max)
    step=$(csg_get_interaction_property step)
    critical rm "${name}.dpot.new"
    do_external table dummy "${min}:${step}:${max}" "${name}.dpot.new"
  fi
}
export -f do_potential_zero_overwrite
for_all "non-bonded" do_potential_zero_overwrite
