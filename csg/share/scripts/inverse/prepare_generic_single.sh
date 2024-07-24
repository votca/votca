#! /usr/bin/env bash
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

show_help () {
  cat <<EOF
${0##*/}, version %version%
This script prepares the potential in step 0, using pot.in or by resampling and
inverting the target distribution

Use --use-table or --use-bi to enforce the method. Otherwise it will use
.pot.in if present and BI if not.

Usage: ${0##*/} [--help] [--use-table|--use-bi]"
EOF
}

use_table=false
use_bi=false
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    --use-table)
        use_table=true
        shift  # past argument
        ;;
    --use-bi)
        use_bi=true
        shift  # past argument
        ;;
    --help)
        show_help
        exit 0
        ;;
    *)
        die "unknown argument $key"
        ;;
    esac
done

if [[ ${use_table} == true && ${use_bi} == true ]]; then
    die "use either --use-table or --use-bi"
fi

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step=$(csg_get_interaction_property step )
comment="$(get_table_comment)"
main_dir=$(get_main_dir)
bondtype="$(csg_get_interaction_property bondtype)"
output="${name}.pot.new"

function table_init() {
    msg "Using given table ${name}.pot.in for ${name}"
    smooth="$(critical mktemp ${name}.pot.in.smooth.XXX)"
    echo "Converting ${main_dir}/${name}.pot.in to ${output}"
    critical csg_resample --in "${main_dir}/${name}.pot.in" --out ${smooth} --grid ${min}:${step}:${max} --comment "$comment"
    extrapolate="$(critical mktemp ${name}.pot.in.extrapolate.XXX)"
    do_external potential extrapolate --type "$bondtype" "${smooth}" "${extrapolate}"
    shifted="$(critical mktemp ${name}.pot.in.shifted.XXX)"
    do_external potential shift --type "${bondtype}" ${extrapolate} ${shifted}
    do_external table change_flag "${shifted}" "${output}"
}

function bi_init() {
    target=$(csg_get_interaction_property inverse.target)
    msg "Using initial guess from dist ${target} for ${name}"
    #resample target dist
    do_external resample target "$(csg_get_interaction_property inverse.target)" "${name}.dist.tgt" 
    # initial guess from rdf
    raw="$(critical mktemp ${name}.pot.new.raw.XXX)"
    kbt="$(csg_get_property cg.inverse.kBT)"
    dist_min="$(csg_get_property cg.inverse.dist_min)"
    do_external dist invert --type "${bondtype}" --kbT "${kbt}" --min "${dist_min}" ${name}.dist.tgt ${raw}
    smooth="$(critical mktemp ${name}.pot.new.smooth.XXX)"
    critical csg_resample --in ${raw} --out ${smooth} --grid ${min}:${step}:${max} --comment "${comment}"
    extrapolate="$(critical mktemp ${name}.pot.new.extrapolate.XXX)"
    do_external potential extrapolate --type "$bondtype" "${smooth}" "${extrapolate}"
    shifted="$(critical mktemp ${name}.pot.new.shifted.XXX)"
    do_external potential shift --type "${bondtype}" ${extrapolate} ${shifted}
    do_external table change_flag "${shifted}" "${output}"
}

table_present=false
if [[ -f ${main_dir}/${name}.pot.in ]]; then
    table_present=true
fi

if [[ $use_bi == true ]]; then
    if [[ $table_present == true ]]; then
        msg "there is a table ${name}.pot.in present, but you still choose BI"
    fi
    bi_init
elif [[ $use_table == true ]]; then
    if [[ $table_present != true ]]; then
        die "missing table ${main_dir}/${name}.pot.in"
    fi
    table_init
else
    # this is the old default behaviour
    if [[ $table_present == true ]]; then
        msg "there is a table ${name}.pot.in present, gonna use it"
        table_init
    else
        bi_init
    fi
fi
