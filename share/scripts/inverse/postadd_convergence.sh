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
postadd convergence script, calcs int of (\${name}.DIST.tgt-\${name}.DIST.new)**2
and saves it to \${name}.conv.
DIST stands for 'dist', but can be changed by onvergence.what option

usage: ${0##*/} infile outfile
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

do_external postadd dummy "$1" "$2"
if [[ $(csg_get_property cg.inverse.method) = "simplex" ]]; then
  msg "WARNING: postadd convergency make no sense for simplex as convergency is calculated anyway - skipping"
  exit 0
fi

name=$(csg_get_interaction_property name)
max=$(csg_get_interaction_property max)
min=$(csg_get_interaction_property min)
step=$(csg_get_interaction_property step)



#these two are arrays
weights=( $(csg_get_interaction_property inverse.post_add_options.convergence.weight 1) )
what_to_do_list=( $(csg_get_interaction_property inverse.post_add_options.convergence.what "dist") )

[[ ${#weights[@]} -ne ${#what_to_do_list[@]} ]] && die "${0##*/}: number of weights does not match number of 'what' to calc convergence from"

sum=0
#we allow multiple thing per interaction to be checked
for ((i=0;i<${#what_to_do_list[@]};i++)); do
  dist=${what_to_do_list[$i]}
  weight=${weights[$i]}
  tmp1="$(critical mktemp ${name}.${dist}.new.resample.XXX)"

  if [[ ! -f "${name}.${dist}.tgt" ]]; then
    #if we need $name.dist.tgt we know how to create it
    if [[ ${dist} = "dist" ]]; then
      do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt"
    else
      die "${0##*/}: file '${name}.${dist}.tgt' was not found. Add the script to create this file to the postadd routine of interaction $name or put it in the maindir and add it to cg.inverse.filelist."
    fi
  fi

  [[ -f ${name}.${dist}.new ]] || die "${0##*/}: file '${name}.${dist}.new' was not found, add a postadd routine of interaction $name to calculate it."
  #resample this, as density dist maybe has the wrong grid
  critical csg_resample --in ${name}.${dist}.new --out $tmp1 --grid "$min:$step:$max"

  diff=$(do_external table combine --sum --no-flags --op d "$tmp1" "${name}.${dist}.tgt")
  is_num "$diff" || die "${0##*/}: strange - result of do_external table difference '$tmp1' and '${name}.${dist}.tgt' was not a number" 
  wdiff=$(csg_calc "$weight" "*" "${diff}")
  echo "Convergence of $dist for ${name} was ${diff} and has weight $weight, so difference is $wdiff"
  sum=$(csg_calc $sum + $wdiff)
done

echo "$sum" > ${name}.conv
