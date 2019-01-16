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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
This script implements the prepares the potential in step 0, using a table pot.in or -1/β log(g) or using the hypernetted chain approximation.

Usage: ${0##*/}
EOF
    exit 0
fi

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
initial_guess="$(csg_get_property cg.inverse.hncgn.initial_guess)"
step=$(csg_get_interaction_property step )
comment="$(get_table_comment)"
main_dir=$(get_main_dir)
bondtype="$(csg_get_interaction_property bondtype)"
output="${name}.pot.new"

case "$initial_guess" in
"table")
    msg "Using given table ${name}.pot.in for ${name}"
    echo "Converting ${main_dir}/${name}.pot.in to ${output}"
    raw="${main_dir}/${name}.pot.in"
    ;;
"bi")
    target=$(csg_get_interaction_property inverse.target)
    msg "Using initial guess from dist ${target} for ${name}"
    #resample target dist
    do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt"
    # initial guess from rdf
    raw="$(critical mktemp ${name}.pot.new.raw.XXX)"
    kbt="$(csg_get_property cg.inverse.kBT)"
    dist_min="$(csg_get_property cg.inverse.dist_min)"
    do_external dist invert --type "${bondtype}" --kbT "${kbt}" --min "${dist_min}" ${name}.dist.tgt ${raw}
    ;;
"hnc")
    target=$(csg_get_interaction_property inverse.target)
    msg "Using initial guess from dist ${target} for ${name}"
    #resample target dist
    do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt"
    # initial guess from rdf with hnc
    raw="$(critical mktemp ${name}.pot.hnc.raw.XXX)"
    kbt="$(csg_get_property cg.inverse.kBT)"
    density="$(csg_get_property cg.inverse.hncgn.density)"
    cut_off="$(csg_get_property cg.inverse.hncgn.cut_off)"
    do_external dist invert_hnc ${name}.dist.tgt "${kbt}" "${density}" "${cut_off}" ${raw}
    ;;
*)
    die "cg.inverse.hncgn.initial_guess has to be either table, bi or hnc"
    ;;
esac

smooth="$(critical mktemp ${name}.pot.new.smooth.XXX)"
critical csg_resample --in ${raw} --out ${smooth} --grid ${min}:${step}:${max} --comment "${comment}"
extrapolate="$(critical mktemp ${name}.pot.new.extrapolate.XXX)"
do_external potential extrapolate --type "$bondtype" "${smooth}" "${extrapolate}"
shifted="$(critical mktemp ${name}.pot.new.shifted.XXX)"
do_external potential shift --type "${bondtype}" ${extrapolate} ${shifted}
do_external table change_flag "${shifted}" "${output}"

