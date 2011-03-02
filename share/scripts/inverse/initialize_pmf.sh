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

pullgroup0=$(csg_get_property cg.non-bonded.pmf.pullgroup0)
pullgroup1=$(csg_get_property cg.non-bonded.pmf.pullgroup1)
conf_init=$(csg_get_property cg.non-bonded.pmf.conf_init)
mdp_init="start_in.mdp"
min=$(csg_get_property cg.non-bonded.pmf.from)
dt=$(get_from_mdp dt "$mdp_init")
rate=$(csg_get_property cg.non-bonded.pmf.rate)
ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
filelist="$(csg_get_property --allow-empty cg.inverse.filelist)"

traj="traj.${ext}"

conf_in=$(csg_get_property --allow-empty cg.non-bonded.pmf.conf_in)
if [! -z "$conf_in" ]; then
    msg "Initial configuration specified in xml, using $conf_in"
    mkdir step_000
    conf_out=$(csg_get_property cg.non-bonded.pmf.conf_in)
    critical cp $(conf_out) step_000/confout.gro
    grompp="$main_dir/grompp.mdp"
    exit 0
fi

# Generate start_in.mdp
cat grompp.mdp.template | sed 's/^pull.*$//' | uniq > tmp
sed -e "s/@TIMESTEP@/$dt/" \
    -e "s/@STEPS@/$steps/" tmp > ${mdp_init}
rm tmp

cp_from_main_dir $filelist
critical cp_from_main_dir grompp.mdp.template ${mdp_init} ${conf_init}  

# Run grompp to generate tpr, then calculate distance
grompp -n index.ndx -c ${conf_init} -f ${mdp_init}
echo -e "${pullgroup0}\n${pullgroup1}" |  g_dist -f ${conf_init} -n index.ndx -o ${conf_init}.xvg
dist=$(sed '/^[#@]/d' ${conf_init}.xvg | awk '{print $2}')
[ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
echo Found distance $dist
cp ${conf_init} conf.gro

# Prepare grommpp file
steps="$(awk "BEGIN{print int(($min-$dist)/$rate/$dt)}")"
if [ $steps -le 0 ]; then
  steps="${steps#-}"
  rate="-$rate"
fi
((steps++))
msg Doing $steps steps with rate $rate

sed -e "s/@DIST@/$dist/" \
    -e "s/@RATE@/$rate/" \
    -e "s/@TIMESTEP@/$dt/" \
    -e "s/@OUT@/0/" \
    -e "s/@PULL_OUT@/0/" \
    -e "s/@STEPS@/$steps/"  grompp.mdp.template > grompp.mdp

# Run simulation to generate initial setup
grompp -n index.ndx
# Create traj so "run gromacs" does not die
touch "$traj"
do_external run gromacs_pmf

# TODO: cleaup as soon as christopg implemented new stuff
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

# Calculate new distance
echo -e "${pullgroup0}\n${pullgroup1}" | g_dist -n index.ndx -f confout.gro
dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
[ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
echo "New distance is $dist"
