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
method="$(csg_get_property cg.inverse.method)"

if [[ -f ${main_dir}/${name}.pot.in ]]; then
  msg "Using given table ${name}.pot.in for ${name}"
  tmp="$(critical mktemp ${name}.pot.in.smooth.XXX)"
  echo "Converting ${main_dir}/${name}.pot.in to ${name}.pot.new through $tmp"
  critical csg_resample --in "${main_dir}/${name}.pot.in" --out ${tmp} --grid ${min}:${step}:${max} --comment "$comment"
  [[ $(csg_get_interaction_property bondtype) = "bonded" ]] && \
  do_external pot shift_bonded ${tmp} ${name}.pot.new || \
  do_external pot shift_nonbonded ${tmp} ${name}.pot.new
else
  [[ $(csg_get_interaction_property bondtype) = "bonded" ]] && die "${0##*/}: Not implemented yet, implement it or provide ${name}.pot.in!"
  target=$(csg_get_interaction_property inverse.target)
  msg "Using initial guess from dist ${target} for ${name}"
  #resample all target dist in $this_dir
  critical csg_resample --in ${main_dir}/${target} --out ${name}.dist.tgt --grid ${min}:${step}:${max} --comment "${comment}"
  if [[ $method = "tf" ]]; then
    #initial guess from density
    do_external calc thermforce ${name}.dist.tgt ${name}.pot.new
  else
    # initial guess from rdf
    tmp="$(critical mktemp ${name}.pot.new.raw.XXX)"
    do_external rdf pot ${name}.dist.tgt ${tmp}
    tmp2="$(critical mktemp ${name}.pot.new.smooth.XXX)"
    critical csg_resample --in ${tmp} --out ${tmp2} --grid ${min}:${step}:${max} --comment "${comment}"
    do_external pot shift_nonbonded ${tmp2} ${name}.pot.new
  fi
fi

