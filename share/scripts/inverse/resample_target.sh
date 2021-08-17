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
This script resamples distribution to grid spacing of the setting xml file and extrapolates if needed

Usage: ${0##*/} input output
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing argument"
input="$1"
main_dir=$(get_main_dir)
[[ -f ${main_dir}/$input ]] || die "${0##*/}: Could not find input file '$input' in maindir ($main_dir)"
output="$2"

min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step=$(csg_get_interaction_property step )
name=$(csg_get_interaction_property name)
tabtype="$(csg_get_interaction_property bondtype)"

comment="$(get_table_comment)"
smooth="$(critical mktemp ${name}.dist.tgt_smooth.XXXXX)"
critical csg_resample --in ${main_dir}/${input} --out ${smooth} --grid ${min}:${step}:${max} --comment "${comment}"
extra="$(critical mktemp ${name}.dist.tgt_extrapolated.XXXXX)"
if [[ $tabtype = "non-bonded" ]]; then
  extra2="$(critical mktemp ${name}.dist.tgt_extrapolated_left.XXXXX)"
  do_external table extrapolate --function linear --avgpoints 1 --region left "${smooth}" "${extra2}"
  # improve RDF in the core region where it is close to zero and sampling is bad using
  # an extrapolation scheme from the better sampled onset of the RDF
  improve_dist_near_core_target="$(csg_get_interaction_property improve_dist_near_core.target)"
  if [[ $improve_dist_near_core_target == "true" ]]; then
      improve_dist_near_core_function="$(csg_get_interaction_property improve_dist_near_core.function)"
      extra3="$(critical mktemp ${name}.dist.tgt_improved.XXXXX)"
      fit_start_g="$(csg_get_interaction_property improve_dist_near_core.fit_start_g)"
      fit_end_g="$(csg_get_interaction_property improve_dist_near_core.fit_end_g)"
      do_external dist improve_near_core --in="${extra2}" --out="${extra3}" \
          --function="$improve_dist_near_core_function" \
          --gmin="$fit_start_g" --gmax="$fit_end_g" 
  elif [[ $improve_dist_near_core_target == "false" ]]; then
      extra3="$extra2"
  else
      die "${0##*/}: improve_dist_near_core.target is ${improve_dist_near_core_target}. Needs to be 'true' or 'false'"
  fi
  do_external table extrapolate --function constant --avgpoints 1 --region right "${extra3}" "${extra}"
elif [[ $tabtype = bond || $tabtype = angle || $tabtype = dihedral ]]; then
  do_external table extrapolate --function linear --avgpoints 1 --region leftright "${smooth}" "${extra}"
else
  die "${0##*/}: Resample of distribution of type $tabtype is not implemented yet"
fi
do_external dist adjust "${extra}" "${output}"
