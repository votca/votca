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
This script implemtents the function post update with the Inverse Boltzmann
Method. This one is not like other post update methods in that the infile is ignored.

Usage: ${0##*/} infile outfile
EOF
   exit 0
fi

name=$(csg_get_interaction_property name)
bondtype="$(csg_get_interaction_property bondtype)"
step_nr=$(get_current_step_nr)
scheme=( $(csg_get_interaction_property inverse.do_potential) )
scheme_nr=$(( (step_nr - 1 ) % ${#scheme[@]} ))
postibi=( $(csg_get_interaction_property inverse.post_update_options.ibi.do) )
postibi_nr=$(( (step_nr - 1 ) % ${#postibi[@]} ))
multistate="$(csg_get_property cg.inverse.multistate.enabled)"


if [[ ${postibi[$postibi_nr]} = 1 ]]; then

  if [[ "${scheme[$scheme_nr]}" == 1 ]]; then
    msg --color blue "WARNING: the potential ${name} has already been updated.
    This ibi post-update will overwrite the original update and all previous
    post-updates! You might want to set do_potential to 0."
  fi

  echo "Apply ibi post-update for interaction ${name}"
  if [[ $multistate == true ]]; then
    state_names_arr=( $(csg_get_property cg.inverse.multistate.state_names) )
    state_weights_arr=( $(csg_get_property cg.inverse.multistate.state_weights) )
    state_kBTs_arr=( $(csg_get_property cg.inverse.multistate.state_kBTs) )
    for s in "${!state_names_arr[@]}"; do
      state="${state_names_arr[s]}"
      weight="${state_weights_arr[s]}"
      kBT="${state_kBTs_arr[s]}"
      is_num "${kBT}" || die "${0##*/}: cg.inverse.multistate.state_kBTs should be numbers, but found '$kBT'"
      # do ibi update per state
      pushd $state
      do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt"
      do_external update ibi_pot "${name}.dist.tgt" "${name}.dist.new" ../"${name}.pot.cur" "${name}.dpot.pure_ibi" "${kBT}"
      # weight each states ibi update
      do_external table scale "${name}.dpot.pure_ibi" "${name}.dpot.ibi_weighted" "${weight}" "${weight}"
      popd
      # sum together
      if (( s == 0 )); then
        cp "${state}/${name}.dpot.pure_ibi" "${name}.dpot.pure_multistate_ibi"
      else
        do_external table combine --no-flags --op "+" "${name}.dpot.pure_multistate_ibi" "${state}/${name}.dpot.pure_ibi" "${name}.dpot.pure_multistate_ibi"
      fi
    done
    do_external potential shift --type "${bondtype}" "${name}.dpot.pure_multistate_ibi" "$2"
  else
    do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt"
    kBT="$(csg_get_property cg.inverse.kBT)"
    is_num "${kBT}" || die "${0##*/}: cg.inverse.kBT should be a number, but found '$kBT'"
    do_external update ibi_pot "${name}.dist.tgt" "${name}.dist.new" "${name}.pot.cur" "${name}.dpot.pure_ibi" "${kBT}"
    do_external potential shift --type "${bondtype}" "${name}.dpot.pure_ibi" "$2"
  fi
else
  echo "No ibi post-update for interaction ${name}"
  do_external postupd dummy "$1" "$2"
fi
