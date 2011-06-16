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
conf="$(csg_get_property cg.inverse.gromacs.conf "conf.gro")"
name=$(csg_get_interaction_property name)
ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
mdp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"
traj="traj.${ext}"
min=$(csg_get_property cg.non-bonded.pmf.min)
mdp_init="$(csg_get_property --allow-empty cg.non-bonded.pmf.mdp_init)"
mdp_prep="$(csg_get_property cg.non-bonded.pmf.mdp_prep)"

cp_from_main_dir "${mdp_prep}"
mv "${mdp_prep}" grompp.mdp

critical cp_from_main_dir ${conf}

#
#Test distance 1st time:
#
# Run grompp to generate tpr, then calculate distance
do_external run gromacs --prepare-only
# Then calcualte the distance
pullgroup0=$(get_from_mdp pull_group0)
pullgroup1=$(get_from_mdp pull_group1)
echo -e "${pullgroup0}\n${pullgroup1}" | g_dist -n index.ndx -f ${conf}
dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
echo "Initial distance is $dist; the expected ditance is $min"

#if the distance between the 2 molecule is not equal to the minumum distance do: 
if `csg_calc ${dist} == ${min}`=="1" do  
echo "Since the 2 distances are not similar the molecules will pe pulled toghether"
critical cp_from_main_dir "${mdp_initi}"
mv "${mdp_prep}" grompp.mdp

index="$(csg_get_property cg.inverse.gromacs.grompp.index "index.ndx")"
[ -f "$index" ] || die "${0##*/}: grompp index file '$index' not found"

# Run simulation to generate initial setup
grompp -n index.ndx ${mdp_opts}
# Create traj so "run gromacs" does not die
touch "$traj"
do_external run gromacs_pmf

done

