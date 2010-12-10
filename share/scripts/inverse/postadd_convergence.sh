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
postadd convergence script, calcs int of (\${name}.DIST.tgt-\${name}.DIST.new)**2
and save it to \${name}.conv. DIST ist dist, but changed by onvergence.what option

usage: ${0##*/} infile outfile

USES: die check_deps do_external wc sed awk paste mktemp

NEEDS: name

OPTIONAL: inverse.post_add_options.convergence.weight inverse.post_add_options.convergence.what
EOF
   exit 0
fi

check_deps "$0"

[ -z "$2" ] && die "${0##*/}: Missing arguments"

do_external postadd dummy "$1" "$2"

name=$(csg_get_interaction_property name)
step=$(csg_get_interaction_property step)

#these two are arrays
weights=( $(csg_get_interaction_property inverse.post_add_options.convergence.weight 1) )
what_to_do_list=( $(csg_get_interaction_property inverse.post_add_options.convergence.what "dist") )

[ ${#weights[@]} -ne ${#what_to_do_list[@]} ] && die "${0##*/}: number of weights does not match number of 'what' to calc convergence from"
tmp="$(true_or_exit mktemp ${name}.conv.XXX)"

#we allow multiple thing per interaction to be checked
for ((i=0;i<${#what_to_do_list[@]};i++)); do
  dist=${what_to_do_list[$i]}
  weight=${weights[$i]}
  tmp1="$(true_or_exit mktemp ${name}.${dist}.tgt.XXX)"
  tmp2="$(true_or_exit mktemp ${name}.${dist}.new.XXX)"
  tmp3="$(true_or_exit mktemp ${name}.${dist}.cmb.XXX)"

  if [ ! -f "${name}.${dist}.tgt" ]; then
    #if we need $name.dist.tgt we know how to create it
    if [ "${dist}" = "dist" ]; then
      do_external resample target
    else
      die "${0##*/}: file '${name}.${dist}.tgt' was not found, add the script to create this file to the postadd routine of interaction $name"
    fi
  fi

  true_or_exit sed -e '/^#/d' -e 's/nan/0.0/g' ${name}.${dist}.tgt > $tmp1
  true_or_exit sed -e '/^#/d' -e 's/nan/0.0/g' ${name}.${dist}.new > $tmp2

  [ $(sed -n '$=' $tmp1) -eq $(sed -n '$=' $tmp2) ] || \
    die "${0##*/}: linenumber of ${name}.${dist}.tgt differs from ${name}.${dist}.new"

  true_or_exit paste $tmp1 $tmp2 > $tmp3
  run_or_exit awk '{if ($4!=$1){print "differ in line NR";exit 1;}}' $tmp3
  echo "Calc convergence for ${name} with weight $weight"
  true_or_exit awk -v bin=$step -v w=$weight -v dist=$dist '{sum+=($5-$2)**2;}END{print dist,sqrt(sum*bin*w);}' $tmp3 >> $tmp
done

true_or_exit awk '{sum+=$2;}END{print sum;}' $tmp > ${name}.conv
