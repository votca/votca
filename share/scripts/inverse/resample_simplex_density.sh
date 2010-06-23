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
This script resamples the target distribution to the grid spacing necessary for calculations.

Usage: ${0##*/} target_directory

USES:  csg_get_interaction_property csg_get_property get_main_dir run_or_exit get_table_comment csg_resample

NEEDS: name min max step target property
EOF
   exit 0
fi

check_deps "$0"

main_dir=$(get_main_dir)
name=$(csg_get_interaction_property name)

comment="$(get_table_comment)"
min=$(csg_get_interaction_property inverse.simplex.density.min)
max=$(csg_get_interaction_property inverse.simplex.density.max)
step=$(csg_get_interaction_property inverse.simplex.density.step)
target=$(csg_get_interaction_property inverse.simplex.density.target)

run_or_exit csg_resample --in ${main_dir}/${target} --out ${name}.dens.tgt --grid ${min}:${step}:${max} --comment "${comment}"
