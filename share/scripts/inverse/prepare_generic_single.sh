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
This script implemtents the function initialize
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES: do_external csg_get_interaction_property log run_or_exit csg_resample log check_deps msg get_main_dir

NEEDS: name min max step

OPTIONAL: inverse.target
EOF
  exit 0
fi

check_deps "$0"

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step=$(csg_get_interaction_property step )
comment="$(get_table_comment)"
main_dir=$(get_main_dir)

if [ -f "${main_dir}/${name}.pot.in" ]; then
  msg "Using given table ${name}.pot.in for ${name}"
  run_or_exit csg_resample --in "${main_dir}/${name}.pot.in" --out ${name}.pot.tmp --grid ${min}:${step}:${max} --comment "$comment"
  do_external pot shift_nb ${name}.pot.tmp ${name}.pot.new
else
  target=$(csg_get_interaction_property inverse.target)
  msg "Using initial guess from dist ${target} for ${name}"
  #copy+resample all target dist in $this_dir
  run_or_exit csg_resample --in ${main_dir}/${target} --out ${name}.dist.tgt --grid ${min}:${step}:${max} --comment "${comment}"
  # RDF_to_POT.pl just does log g(r) + extrapolation
  do_external rdf pot ${name}.dist.tgt ${name}.pot.new
fi

