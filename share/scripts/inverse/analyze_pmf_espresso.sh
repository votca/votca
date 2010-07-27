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
# Unless required by applicale law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if [ "$1" = "--help" ]; then
    cat <<EOF
${0##*/}, version %version%
This script analyzes the pmf for espresso

Usage: ${0##*/}

USES: csg_get_interaction_property csg_get_property log run_or_exit csg_resample is_done mark_done msg check_deps csg_part_dist

NEEDS: type cg.inverse.espresso.min cg.inverse.espresso.binsize cg.inverse.espresso.max

OPTIONAL: cg.inverse.espresso.current_pmf cg.inverse.espresso.target_pmf cg.inverse.espresso.traj cg.inverse.espresso.analysis
EOF
    exit 0
fi

check_deps "$0"

# Topology+Trajectory
toptraj="$(csg_get_property cg.inverse.espresso.traj "top_traj.esp")"
[ -f "$toptraj" ] || die "${0##*/}: Topology+Trajectory file $toptraj not found"

# Particle distribution file (temp)
part_dist="$(mktemp part.dist.XXXXX)"

# Current PMF
pmf_cur="$(csg_get_property cg.inverse.espresso.current_pmf "current.pmf.dat")"
[ -f "$pmf_cur" ] || die "${0##*/}: Current pmf $pmf_cur not found"

# Target PMF
pmf_tgt="$(csg_get_property cg.inverse.espresso.target_pmf "target.pmf.dat")"
[ -f "$pmf_tgt" ] || die "${0##*/}: Target pmf $pmf_tgt not found"

# particle types involved in the optimization (except inclusion)
ptypes="$(for_all non-bonded csg_get_interaction_property type)"
# write to file
file_ptypes="$(mktemp file.ptypes.XXXXX)"
echo $ptypes > $file_ptypes
num_ptypes="$(cat $file_ptypes | wc -w)"

# Grid for PMF
min="$(csg_get_property cg.inverse.espresso.min)"
[ -z "$min" ] && die "${0##*/}: Could not read espresso property min"
binsize="$(csg_get_property cg.inverse.espresso.binsize)"
[ -z "$min" ] && die "${0##*/}: Could not read espresso property binsize"
max="$(csg_get_property cg.inverse.espresso.max)"
[ -z "$min" ] && die "${0##*/}: Could not read espresso property max"

# Output of pmf analysis
err_ptypes="$(csg_get_property cg.inverse.espresso.analysis "err.ptypes.dat")"

log "Analyzing PMF"
if is_done "pmf-analysis"; then
    msg "PMF analysis is already done"
else

    run_or_exit csg_part_dist --top $toptraj --trj $toptraj --out $part_dist --grid ${min}:${binsize}:${max} --ptypes $file_ptypes --shift_com

    comment="$(get_table_comment)"
    run_or_exit csg_resample --in ${pmf_cur} --out ${pmf_cur}.new --grid ${min}:${binsize}:${max} --comment "$comment"
    run_or_exit csg_resample --in ${pmf_tgt} --out ${pmf_tgt}.new --grid ${min}:${binsize}:${max} --comment "$comment"

    run_or_exit do_external pmf analyze ${pmf_cur}.new ${pmf_tgt}.new $part_dist $num_ptypes $err_ptypes 

    mark_done "pmf-analysis"
fi
