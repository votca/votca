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

#if user gave us a conf_in we use this one
conf_in=$(csg_get_interaction_property --allow-empty pmf.conf_in)
if [[ -n $conf_in ]]; then
    msg "Initial configuration specified in xml, using $conf_in"
    conf_out="$(csg_get_property cg.inverse.gromacs.conf_out "confout.gro")"
    cp_from_main_dir --rename "${conf_in}" "${conf_out}"
    #todo check if distance is really correct
    exit 0
fi

name=$(csg_get_interaction_property name)
ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
mdp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"
traj="traj.${ext}"

# Generate start_in.mdp
mdp_template="$(csg_get_property cg.inverse.gromacs.mdp_template "grompp.mdp.template")"
mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
cp_from_main_dir --rename "${mdp_template}" "${mdp}"
sed -e "s/@STEPS@/$steps/" \
    -e "s/@EXCL@//" \
    -e "s/@OUT@/0/" -i  "${mdp}"

conf=$(csg_get_interaction_property cg.inverse.gromacs.conf "conf.gro")
critical cp_from_main_dir ${conf} 

# Run grompp to generate tpr, then calculate distance
do_external run gromacs --prepare-only

grp1=$(csg_get_interaction_property --allow-empty gromacs.grp1)
grp2=$(csg_get_interaction_property --allow-empty gromacs.grp2)
g_dist="$(csg_get_property cg.inverse.gromacs.g_dist.bin "g_dist")"
[ -n "$(type -p $grompp)" ] || die "${0##*/}: grompp binary '$grompp' not found"
index="$(csg_get_property cg.inverse.gromacs.grompp.index "index.ndx")"
[ -f "$index" ] || die "${0##*/}: grompp index file '$index' not found"

distfile="$(critical mktemp smooth_${name}.dist.XXXXX)"
echo -e "${grp1}\n${grp2}" |  "$g_dist" -f "${conf}" -n "$index" -o "$distfile"
dist=$(awk '/^[^#@]/{print $2}' "$distfile") || die "${0##*/}: Awk get of distance from $distfile failed"
[ -z "$dist" ] && die "${0##*/}: Could not fetch dist for $distfile"
echo "Distance between ${grp1} and ${grp2} in ${conf}:$dist"

#ToDo

# Prepare grommpp file
dt=$(get_simulation_setting dt)
rate=$(csg_get_interaction_property pmf.rate)
min=$(csg_get_interaction_property min)
steps="$(awk "BEGIN{print int(($min-$dist)/$rate/$dt)}")"
if [ $steps -le 0 ]; then
  steps="${steps#-}"
  rate="-$rate"
fi
((steps++))

msg Doing $steps steps with rate $rate

sed -e "s/@DIST@/$dist/" \
    -e "s/@RATE@/$rate/" \
    -e "s/@STEPS@/$steps/" \
    -e "s/@EXCL@//" \
    -e "s/@OUT@/0/"  grompp.mdp.template > grompp.mdp

# Run simulation to generate initial setup
grompp -n index.ndx ${mdp_opts}
# Create traj so "run gromacs" does not die
touch "$traj"
do_external run gromacs_pmf

# TODO: cleaup as soon as christopg implemented new stuff
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

# Calculate new distance
echo -e "${pullgroup0}\n${pullgroup1}" | g_dist -n index.ndx -f confout.gro
dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
[ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
echo "New distance is $dist"
