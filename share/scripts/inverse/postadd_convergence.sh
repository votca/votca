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
DIST is dist, but changed by onvergence.what option

usage: ${0##*/} infile outfile
EOF
   exit 0
fi

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"

do_external postadd dummy "$1" "$2"

name=$(csg_get_interaction_property name)
max=$(csg_get_interaction_property max)
min=$(csg_get_interaction_property min)
step=$(csg_get_interaction_property step)

#these two are arrays
weights=( $(csg_get_interaction_property inverse.post_add_options.convergence.weight 1) )
what_to_do_list=( $(csg_get_interaction_property inverse.post_add_options.convergence.what "dist") )

[ ${#weights[@]} -ne ${#what_to_do_list[@]} ] && die "${0##*/}: number of weights does not match number of 'what' to calc convergence from"
tmp="$(critical mktemp ${name}.conv.XXX)"

#we allow multiple thing per interaction to be checked
for ((i=0;i<${#what_to_do_list[@]};i++)); do
  dist=${what_to_do_list[$i]}
  weight=${weights[$i]}
  tmp1="$(critical mktemp ${name}.${dist}.tgt.XXX)"
  tmp2="$(critical mktemp ${name}.${dist}.newcut.XXX)"
  tmp3="$(critical mktemp ${name}.${dist}.new.XXX)"
  tmp4="$(critical mktemp ${name}.${dist}.cmb.XXX)"

  if [ ! -f "${name}.${dist}.tgt" ]; then
    #if we need $name.dist.tgt we know how to create it
    if [ "${dist}" = "dist" ]; then
      do_external resample target
    else
      die "${0##*/}: file '${name}.${dist}.tgt' was not found, add the script to create this file to the postadd routine of interaction $name"
    fi
  fi

  critical sed -e '/^#/d' -e 's/nan/0.0/g' ${name}.${dist}.tgt > $tmp1
  critical csg_resample --in ${name}.${dist}.new --out $tmp2 --grid "$min:$step:$max"
  critical sed -e '/^#/d' -e 's/nan/0.0/g' $tmp2 > $tmp3

  [ $(sed -n '$=' $tmp1) -eq $(sed -n '$=' $tmp3) ] || \
    die "${0##*/}: linenumber of ${name}.${dist}.tgt differs from ${name}.${dist}.new"

  critical paste $tmp1 $tmp3 > $tmp4
  critical awk '{if ($4!=$1){print "x column differs in line",NR;exit 1;}}' $tmp4
  echo "Calc convergence for ${name} with weight $weight"
  critical awk -v bin=$step -v w=$weight -v dist=$dist '{sum+=($5-$2)**2;}END{print dist,sqrt(sum*bin*w);}' $tmp4 >> $tmp
done

critical awk '{sum+=$2;}END{print sum;}' $tmp > ${name}.conv
