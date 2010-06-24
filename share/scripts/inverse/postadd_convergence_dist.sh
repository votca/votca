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
postadd convergence script, calcs int of (\${name}.dist.tgt-\${name}.dist.new)**2
and save it to \${name}.conv

Usage: ${0##*/} infile outfile

USES: die check_deps do_external wc sed awk paste mktemp

NEEDS: name

OPTIONAL: inverse.post_add_options.convergence.weight
EOF
   exit 0
fi

check_deps "$0"

[ -z "$2" ] && die "${0##*/}: Missing arguments"

do_external postadd dummy "$1" "$2"

name=$(csg_get_interaction_property name)
step=$(csg_get_interaction_property step)
weight=$(csg_get_interaction_property inverse.post_add_options.convergence.weight 1)

tmp1="$(true_or_exit mktemp ${name}.dist.tgt.XXX)"
tmp2="$(true_or_exit mktemp ${name}.dist.new.XXX)"
tmp3="$(true_or_exit mktemp ${name}.dist.cmb.XXX)"

if [ ! -f "${name}.dist.tgt" ]; then
  do_external resample target
fi

true_or_exit sed -e '/^#/d' -e 's/nan/0.0/g' ${name}.dist.tgt > $tmp1
true_or_exit sed -e '/^#/d' -e 's/nan/0.0/g' ${name}.dist.new > $tmp2

[ $(sed -n '$=' $tmp1) -eq $(sed -n '$=' $tmp2) ] || \
  die "${0##*/}: linenumber of ${name}.dist.tgt differs from ${name}.dist.new"

true_or_exit paste $tmp1 $tmp2 > $tmp3
run_or_exit awk '{if ($4!=$1){print "differ in line NR";exit 1;}}' $tmp3
log "Calc convergence for ${name} with weight $weight"
true_or_exit awk -v bin=$step -v w=$weight '{sum+=($5-$2)**2;}END{print sqrt(sum*bin*w);}' $tmp3 > ${name}.conv

