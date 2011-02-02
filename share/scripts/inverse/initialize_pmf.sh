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

USES: do_external csg_get_interaction_property check_deps 

NEEDS: pullgroup0 pullgroup1 pullgroup0_type pullgroup1_type confin min max step dt rate kB
EOF
  exit 0
fi

check_deps "$0"

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

# Get atom numbers for pull_group0, pull_group1
n1="$(get_pullgroup_nr "$pullgroup0" "$conf_in.gro")"
n2="$(get_pullgroup_nr "$pullgroup1" "$conf_in.gro")"
if [ $n1 -eq $n2 ]; then
  n2="$(get_pullgroup_nr "$pullgroup1" "$conf_in.gro" 2)"
fi

# Generate index files for pullgroup0, pullgroup1 and their environments
echo -e "
del 1-10
r SOL
a $n1
a $n2
name 2 pullgroup0
name 3 pullgroup1
a $pullgroup0_type & ! a $n1 & ! a $n2
name 4 $pullgroup0_type
a $pullgroup1_type & ! a $n1 & ! a $n2
name 5 $pullgroup1_type
q" | run make_ndx -f ${conf_in}.gro

# Run grompp to generate tpr, then calculate distance
run grompp -n index.ndx -c ${conf_in}.gro -o ${conf_in}.tpr -f ${conf_in}.mdp -po ${conf_in}_all.mdp
echo -e "pullgroup0\npullgroup1" | run g_dist -f ${conf_in}.gro -s ${conf_in}.tpr -n index.ndx -o ${conf_in}.xvg
dist=$(sed '/^[#@]/d' ${conf_in}.xvg | awk '{print $2}')
[ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
echo Found distance $dist
cp $confin.gro conf.gro

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
    -e "s/@STEPS@/$steps/" grompp.mdp.template > grompp.mdp

# Run simulation to generate initial setup
run --log log_grompp2 grompp -n index.ndx
do_external run gromacs_pmf

# Calculate new distance
echo -e "pullgroup0\npullgroup1" | run g_dist -n index.ndx -f confout.gro
dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
[ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
echo "New distance is $dist"

echo "Initial setup is done"
exit 0
