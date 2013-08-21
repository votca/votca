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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script implements the prepares the potential in step 0, using pot.in or by resampling the target distribution

Usage: ${0##*/}
EOF
  exit 0
fi

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step=$(csg_get_interaction_property step )
comment="$(get_table_comment)"
main_dir=$(get_main_dir)
bondtype="$(csg_get_interaction_property bondtype)"
output="${name}.pot.new"

if [[ -f ${main_dir}/${name}.pot.in ]]; then
  msg "Using given table ${name}.pot.in for ${name}"
  smooth="$(critical mktemp ${name}.pot.in.smooth.XXX)"
  echo "Converting ${main_dir}/${name}.pot.in to ${output}"
  critical csg_resample --in "${main_dir}/${name}.pot.in" --out ${smooth} --grid ${min}:${step}:${max} --comment "$comment"
  extrapolate="$(critical mktemp ${name}.pot.in.extrapolate.XXX)"
  do_external potential extrapolate --type "$bondtype" "${smooth}" "${extrapolate}"
  do_external table change_flag "${extrapolate}" "${output}"
else
  target=$(csg_get_interaction_property inverse.target)
  msg "Using initial guess from dist ${target} for ${name}"
  if [[ $bondtype = "thermforce" ]]; then
    #therm force is resampled later and as one want to symetrize 1d density
    cp_from_main_dir --rename "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt" 
    #initial guess from density
    raw="$(critical mktemp -u ${name}.pot.new.raw.XXX)"
    do_external calc thermforce ${name}.dist.tgt ${raw}
    do_external table change_flag "${raw}" "${output}"
  elif [[ ${bondtype} = "non-bonded" ]]; then
    #resample target dist
    do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt" 
    # initial guess from rdf
    raw="$(critical mktemp ${name}.pot.new.raw.XXX)"
    kbt="$(csg_get_property cg.inverse.kBT)"
    dist_min="$(csg_get_property cg.inverse.rdf_min)"
    do_external dist invert --type "${bondtype}" --kbT "${kbt}" --min "${dist_min}" ${name}.dist.tgt ${raw}
    smooth="$(critical mktemp ${name}.pot.new.smooth.XXX)"
    critical csg_resample --in ${raw} --out ${smooth} --grid ${min}:${step}:${max} --comment "${comment}"
    extrapolate="$(critical mktemp ${name}.pot.new.extrapolate.XXX)"
    do_external potential shift --type "${bondtype}" ${smooth} ${extrapolate}
    do_external table change_flag "${extrapolate}" "${output}"
  else
    die "${0##*/}: Not implemented yet, implement it or provide ${name}.pot.in!"
  fi
fi

