#! /bin/bash -e
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
This script implemtents the function initialize for the PMF calculator

Usage: ${0##*/}

EOF
  exit 0
fi

mdp_init=$(csg_get_property --allow-empty cg.non-bonded.pmf.mdp_init)
mdp_prep=$(csg_get_property --allow-empty cg.non-bonded.pmf.mdp_prep)

# Get grompp file
if [[ -z "$mdp_init" ]]; then  
  cp_from_main_dir conf.gro
  cp_from_main_dir ${mdp_prep}
  mv ${mdp_prep} grompp.mdp
else
  cp_from_main_dir conf.gro
  cp_from_main_dir ${mdp_init}
  mv ${mdp_init} grompp.mdp
fi

pullgroup0=$(get_simulation_setting pull_group0)
pullgroup1=$(get_simulation_setting pull_group1)
steps=$(get_simulation_setting nsteps)
rate=$(get_simulation_setting pull_rate1)
dt=$(get_simulation_setting dt)

min=$(csg_get_property cg.non-bonded.pmf.min)
mdp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"
index=$(csg_get_property cg.inverse.gromacs.grompp.index "index.ndx")
[ -f "$index" ] || die "${0##*/}: grompp index file '$index' not found"
ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
traj="traj.${ext}"

# Check if distance between pullgroups is at min
do_external run gromacs --prepare-only
echo -e "${pullgroup0}\n${pullgroup1}" | g_dist -n ${index} -f conf.gro
dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
msg "Initial distance between pullgroups is $dist"

if [[ `csg_calc ${dist} == ${min}`=="1" ]]; then
  msg "Now pulling to distance $min"
  [ -n "$mdp_init" ] || die "${0##*/}: initialization grompp file is not provided"
  steps="$(awk "BEGIN{print int(($min-$dist)/$rate/$dt)}")"
  # Determine whether to pull apart or together
  if [ $steps -le 0 ]; then
    steps="${steps#-}"
    rate="-$rate"
  fi
  ((steps++))
  sed -i -e "s/nsteps.*$/nsteps                   = $steps/" \
         -e "s/pull_rate1.*$/pull_rate1               = $rate/" \
         -e "s/pull_init1.*$/pull_init1               = $dist/" grompp.mdp
  grompp -n ${index} ${mdp_opts}
  do_external run gromacs_pmf
else
  msg "Pullgrous are already at $min distance"
fi
