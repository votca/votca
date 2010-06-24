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
This script:
- reads simplex infile and creates a table (prepare simplex)
- calculates potential for first parameter set (pot [function])
- flags that parameter set as 'active' (par pot)

Usage: ${0##*/}

USES: do_external csg_get_interaction_property cp_from_main_dir run_or_exit check_deps msg

NEEDS: name
EOF
  exit 0
fi

check_deps "$0"

main_dir=$(get_main_dir)
name=$(csg_get_interaction_property name)
property=$(csg_get_property cg.inverse.simplex.property)
count=$(echo "$property" | wc -w)
function=$(csg_get_interaction_property inverse.simplex.function)
param_N=$(do_external pot $function --nparams)

if [ -f $main_dir/simplex_$name.in ]; then
  msg "Converting table simplex_${name}.in for ${name}"
  cp_from_main_dir simplex_$name.in
  do_external prepare_table simplex simplex_$name.in simplex_$name.cur state_$name.new $param_N
else
  die "No input file simplex_$name.in found"
fi

cp simplex_$name.cur simplex_$name.tmp
if [ $count -gt "1" ]; then
  for p in $property; do
    cp simplex_$name.tmp simplex_$name\_$p.tmp
  done
fi

msg "Calculating potential from first initial guess"
run_or_exit do_external par pot simplex_$name.cur simplex_$name.new $param_N 0
