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
This script implements update step of relative entropy method by csg_reupdate program

Usage: ${0##*/}
EOF
   exit 0
fi

topol=$(csg_get_property --allow-empty cg.inverse.gromacs.re.topol)
[[ -z ${topol} ]] && topol="$(csg_get_property cg.inverse.gromacs.topol_out)"
[[ -f $topol ]] || die "${0##*/}: gromacs topol file '$topol' not found, possibly you have to add it to cg.inverse.filelist" 
ext=$(csg_get_property cg.inverse.gromacs.traj_type)
traj="traj.${ext}"
[[ -f $traj ]] || die "${0##*/}: gromacs traj file '$traj' not found"

csg_reupdate_opts="$(csg_get_property --allow-empty cg.inverse.re.csg_reupdate.opts)"

if is_done "re_update"; then
  echo "re update is already done"
else
  critical csg_reupdate --top ${topol} --trj $traj --options $CSGXMLFILE ${csg_reupdate_opts}
  mark_done "re_update"
fi

if [[ -f "notsympos" ]]; then
  msg "re updated using the steepest descent"
fi
