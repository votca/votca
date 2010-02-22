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
This script resamples target distribution to grid spacing
for calculations

Usage: ${0##*/} target_directory

USES:  csg_get_interaction_property run_or_exit csg_resample check_deps get_main_dir

NEEDS: min max step inverse.target name
EOF
   exit 0
fi

check_deps "$0"

min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step=$(csg_get_interaction_property step )
target=$(csg_get_interaction_property inverse.target)
name=$(csg_get_interaction_property name)
main_dir=$(get_main_dir)

run_or_exit csg_resample --in ${main_dir}/${target} --out ${name}.dist.tgt --grid ${min}:${step}:${max}
