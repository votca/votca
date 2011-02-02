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

if [ -f "${main_dir}/${name}.pot.in" ]; then
  msg "Using given table ${name}.pot.in for ${name}"
  tmp="$(critical mktemp ${name}.pot.in.smooth.XXX)"
  echo "Converting ${main_dir}/${name}.pot.in to ${name}.pot.new through $tmp"
  critical csg_resample --in "${main_dir}/${name}.pot.in" --out ${tmp} --grid ${min}:${step}:${max} --comment "$comment"
  do_external pot shift_nb ${tmp} ${name}.pot.new
else
  target=$(csg_get_interaction_property inverse.target)
  msg "Using initial guess from dist ${target} for ${name}"
  #copy+resample all target dist in $this_dir
  critical csg_resample --in ${main_dir}/${target} --out ${name}.dist.tgt --grid ${min}:${step}:${max} --comment "${comment}"
  # RDF_to_POT.pl just does log g(r) + extrapolation
  do_external rdf pot ${name}.dist.tgt ${name}.pot.new
fi

