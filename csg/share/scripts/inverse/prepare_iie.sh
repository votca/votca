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
This script calculates dc/dh once for all iie steps

Usage: ${0##*/}
EOF
   exit 0
fi

do_external prepare generic

tgt_dcdh="$(csg_get_property cg.inverse.iie.tgt_dcdh)"
if [[ $tgt_dcdh == 'true' ]]; then
  if [[ -f ${get_main_dir}/dcdh.npz ]]; then
    msg "dcdh.npz is already present in main dir, using it."
    exit 0
  fi
  msg "Calculating dc/dh for all later iterations"
  # make sure dist and dist-intra are here
  for_all "non-bonded" do_external resample target '$(csg_get_interaction_property inverse.target)' '$(csg_get_interaction_property name).dist.tgt'
  if [[ $(csg_get_property cg.inverse.initial_guess.method) != ie ]]; then
    # resample intramolecular only if present. Later iie.py will only load the ones that are needed
    for_all "non-bonded" do_external resample target --no-extrap --skip-if-missing '$(csg_get_interaction_property inverse.target_intra)' '$(csg_get_interaction_property name).dist-intra.tgt'
  fi
  # verbose
  verbose=$(csg_get_property cg.inverse.initial_guess.ie.verbose)
  step_nr=$(get_current_step_nr)
  [[ "${verbose}" == 'true' ]] && verbose_flag="--verbose"
  [[ "${verbose}" == 'step0+1' ]] && [[ $step_nr == '0' || $step_nr == '1' ]] && verbose_flag="--verbose"
  # topology for molecular conections and volume
  topol=$(csg_get_property cg.inverse.iie.topol)
  [[ -f $topol ]] || die "${0##*/}: topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist"
  # volume
  volume=$(critical csg_dump --top "$topol" | grep 'Volume' | awk '{print $2}')
  ([[ -n "$volume" ]] && is_num "$volume") || die "could not determine the volume from file ${topol}"
  do_external dist invert_iie dcdh \
  ${verbose_flag-} \
  --volume "$volume" \
  --topol "$topol" \
  --options "$CSGXMLFILE" \
  --g-tgt-ext "dist.tgt" \
  --g-tgt-intra-ext "dist-intra.tgt" \
  --out ${get_main_dir}/dcdh.npz
fi
