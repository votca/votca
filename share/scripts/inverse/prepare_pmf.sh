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
This script implemtents the function prepare for the PMF calculator

Usage: ${0##*/}

EOF
  exit 0
fi

conf_start="start"
pullgroup0=$(csg_get_property cg.non-bonded.pmf.pullgroup0)
pullgroup1=$(csg_get_property cg.non-bonded.pmf.pullgroup1)
min=$(csg_get_property cg.non-bonded.pmf.from)
max=$(csg_get_property cg.non-bonded.pmf.to)
rate=$(csg_get_property cg.non-bonded.pmf.rate)
mdp_init="start_in.mdp"
mdp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"

# Generate start_in.mdp
critical cp_from_main_dir grompp.mdp.template
cat grompp.mdp.template | sed 's/^pull.*$//' | uniq > tmp
sed -e "s/@STEPS@/$steps/" \
    -e "s/@EXCL@//" \
    -e "s/@OUT@/1/" tmp > ${mdp_init}
rm tmp

dt=$(get_from_mdp dt "$mdp_init")

# Run grompp to generate tpr, then calculate distance
grompp -n index.ndx -c conf.gro -o ${conf_start}.tpr -f start_in.mdp -po ${conf_start}.mdp ${mdp_opts}
# TODO: do not hardcode pullgroups
echo -e "${pullgroup0}\n${pullgroup1}" | g_dist -f conf.gro -s ${conf_start}.tpr -n index.ndx -o ${conf_start}.xvg

# TODO: check for all commands (sed, awk, ...) whether successful or use bash -e
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
msg Doing $(($steps+1)) simulations with rate $rate, output every $out steps at $dt ps 
sed -e "s/@DIST@/$dist/" \
    -e "s/@RATE@/$rate/" \
    -e "s/@STEPS@/$steps/" \
    -e "s/@EXCL@//" \
    -e "s/@OUT@/1/" grompp.mdp.template > grompp.mdp

# Run simulation to generate initial setup
grompp -n index.ndx ${mdp_opts}
do_external run gromacs_pmf

# Wait for job to finish when running in background
confout="$(csg_get_property cg.inverse.gromacs.conf_out "confout.gro")"
background=$(csg_get_property --allow-empty cg.inverse.simulation.background "no")
sleep_time=$(csg_get_property --allow-empty cg.inverse.simulation.sleep_time "60")
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
echo -e "${pullgroup0}\n${pullgroup1}" | g_dist -n index.ndx
echo "System" | trjconv -sep -o conf_start.gro 
