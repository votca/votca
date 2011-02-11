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
This script implemtents the function prepare for the PMF calculator

Usage: ${0##*/}

EOF
  exit 0
fi

if [ -f ${CSGSHARE}/scripts/inverse/functions_pmf.sh ]; then
  source ${CSGSHARE}/scripts/inverse/functions_pmf.sh || die "Could not source functions_pmf.sh"
else
  die "Could not find functions_pmf.sh"
fi

pullgroup0=$(csg_get_property cg.non-bonded.name | sed 's/-.*$//')
pullgroup1=$(csg_get_property cg.non-bonded.name | sed 's/^.*-//')
conf_start="start"
min=$(csg_get_property cg.non-bonded.min)
max=$(csg_get_property cg.non-bonded.max)
dt=$(csg_get_property cg.non-bonded.dt)
rate=$(csg_get_property cg.non-bonded.rate)

# Run grompp to generate tpr, then calculate distance
run grompp -n index.ndx -c conf.gro -o ${conf_start}.tpr -f start_in.mdp -po ${conf_start}.mdp
echo -e "pullgroup0\npullgroup1" | run g_dist -f conf.gro -s ${conf_start}.tpr -n index.ndx -o ${conf_start}.xvg
dist=$(sed '/^[#@]/d' ${conf_start}.xvg | awk '{print $2}')
[ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
echo Found distance $dist

# Prepare grompp file
steps="$(awk "BEGIN{print int(($max-$dist)/$rate/$dt)}")"
if [ $steps -le 0 ]; then
  steps="${steps#-}"
  rate="-$rate"
fi
((steps++))
out="$(awk "BEGIN{print int(1/$dt)}")"
echo Doing $steps steps with rate $rate output every $out steps a $dt ps 
sed -e "s/@DIST@/$dist/" \
    -e "s/@RATE@/$rate/" \
    -e "s/@TIMESTEP@/$dt/" \
    -e "s/@OUT@/$out/" \
    -e "s/@PULL_OUT@/0/" \
    -e "s/@STEPS@/$steps/" grompp.mdp.template > grompp.mdp

# Run simulation to generate initial setup
run --log log_grompp2 grompp -n index.ndx
do_external run gromacs_pmf

# Wait for job to finish when running in background
confout="$(csg_get_property cg.inverse.gromacs.conf_out "confout.gro")"
background=$(csg_get_property --allow-empty cg.inverse.parallel.background "no")
sleep_time=$(csg_get_property --allow-empty cg.inverse.parallel.sleep_time "60")
sleep 10
if [ "$background" == "yes" ]; then
  while [ ! -f "$confout" ]; do
    sleep $sleep_time
  done
else
  ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
  traj="traj.${ext}"
  [ -f "$confout" ] || die "${0##*/}: Gromacs end coordinate '$confout' not found after running mdrun"
fi

# Calculate new distance and divide trj into separate frames
echo -e "pullgroup0\npullgroup1" | run g_dist -n index.ndx
echo "System" | trjconv -sep -o conf_start.gro 
