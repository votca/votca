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
if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script resamples target distribution to grid spacing of the setting xml file

Usage: ${0##*/}
EOF
   exit 0
fi

min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step=$(csg_get_interaction_property step )
target=$(csg_get_interaction_property inverse.target)
name=$(csg_get_interaction_property name)
main_dir=$(get_main_dir)
output="${name}.dist.tgt"

comment="$(get_table_comment)"
smooth="$(critical mktemp ${name}.dist.tgt_smooth.XXXXX)"
critical csg_resample --in ${main_dir}/${target} --out ${smooth} --grid ${min}:${step}:${max} --comment "${comment}"

tabtype="$(csg_get_interaction_property bondtype)"
extrapol="$(critical mktemp ${name}.dist.tgt_extrapol.XXXXX)"
if [[ $tabtype = "non-bonded" || $tabtype = "C6" || $tabtype = "C12" ]]; then
  #the left side is usually not a problem, but still we do it
  do_external table extrapolate --function constant --avgpoints 1 --region leftright "${smooth}" "${output}"
else
  die "${0##*/}: Resample of bonded distribution is not implemented yet"
  #exponential can lead to values smaller 0, needs to be checked again
  do_external table extrapolate --function exponential --avgpoints 1 --region leftright "${smooth}" "${output}"
fi
