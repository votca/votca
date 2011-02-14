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
This script implemtents the function initialize for the PMF calculator

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
pullgroup0_type=$(csg_get_property cg.non-bonded.type1)
pullgroup1_type=$(csg_get_property cg.non-bonded.type2)
conf_in=$(csg_get_property cg.non-bonded.conf_in)
min=$(csg_get_property cg.non-bonded.min)
max=$(csg_get_property cg.non-bonded.max)
dt=$(csg_get_property cg.non-bonded.dt)
rate=$(csg_get_property cg.non-bonded.rate)
ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
traj="traj.${ext}"

# Get atom numbers for pull_group0, pull_group1
n1="$(get_pullgroup_nr "$pullgroup0" "$conf_in.gro")"
n2="$(get_pullgroup_nr "$pullgroup1" "$conf_in.gro")"
if [ $n1 -eq $n2 ]; then
  n2="$(get_pullgroup_nr "$pullgroup1" "$conf_in.gro" 2)"
fi

# Generate index files for pullgroup0, pullgroup1 and their environments
echo -e "
del 1-10
a $n1
a $n2
name 1 pullgroup0
name 2 pullgroup1
a $pullgroup0_type & ! a $n1 & ! a $n2
name 3 pullgroup0_type
a $pullgroup1_type & ! a $n1 & ! a $n2
name 4 pullgroup1_type
q" | run make_ndx -f ${conf_in}.gro

# Get remaining energy groups
last=$(wc -l ${conf_in}.gro | awk '{print $1}')
sec_last=$(($last-1))
energygrps="$(sed -n "3,${sec_last}p" ${conf_in}.gro | awk '{print $1}' | sed 's/[0-9]//g' | sort | uniq | grep -v "$pullgroup0_type" | grep -v "$pullgroup1_type" | xargs)"

for i in $energygrps; do
  echo -e "
r $i
q" | run make_ndx -f ${conf_in}.gro -n index.ndx -o index2.ndx
done

# Run grompp to generate tpr, then calculate distance
run grompp -n index.ndx -c ${conf_in}.gro -o ${conf_in}.tpr -f ${conf_in}.mdp -po ${conf_in}_all.mdp
echo -e "pullgroup0\npullgroup1" | run g_dist -f ${conf_in}.gro -s ${conf_in}.tpr -n index.ndx -o ${conf_in}.xvg
dist=$(sed '/^[#@]/d' ${conf_in}.xvg | awk '{print $2}')
[ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
echo Found distance $dist
cp ${conf_in}.gro conf.gro

# Prepare grommpp file
steps="$(awk "BEGIN{print int(($min-$dist)/$rate/$dt)}")"
if [ $steps -le 0 ]; then
  steps="${steps#-}"
  rate="-$rate"
fi
((steps++))
echo Doing $steps steps with rate $rate
sed -e "s/@DIST@/$dist/" \
    -e "s/@RATE@/$rate/" \
    -e "s/@TIMESTEP@/$dt/" \
    -e "s/@OUT@/0/" \
    -e "s/@PULL_OUT@/0/" \
    -e "s/@STEPS@/$steps/" \
    -e "s/@ENERGYGRPS/pullgroup0 pullgroup1 pullgroup0_type pullgroup1_type $energygrps" grompp.mdp.template > grompp.mdp

# Run simulation to generate initial setup
run --log log_grompp2 grompp -n index.ndx
# Create traj so "run gromacs" does not die
touch "$traj"
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

# Calculate new distance
echo -e "pullgroup0\npullgroup1" | run g_dist -n index.ndx -f confout.gro
dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
[ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
echo "New distance is $dist"
