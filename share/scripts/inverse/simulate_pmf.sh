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

USES: do_external csg_get_interaction_property check_deps

NEEDS: pullgroup0 pullgroup1 confin min max step dt rate kB
EOF
  exit 0
fi

# Read previous grompp file to define pullgroups
mdp_prep=$(csg_get_property cg.non-bonded.pmf.mdp_prep)
cp_from_main_dir $mdp_prep
mv $mdp_prep grompp.mdp

pullgroup0=$(get_simulation_setting pull_group0)
pullgroup1=$(get_simulation_setting pull_group1)

# Overwrite with new grompp file
mdp_sim=$(csg_get_property cg.non-bonded.pmf.mdp_sim)
cp_from_main_dir $mdp_sim
mv $mdp_sim grompp.mdp

min=$(csg_get_property cg.non-bonded.pmf.min)
max=$(csg_get_property cg.non-bonded.pmf.max)
filelist="$(csg_get_property --allow-empty cg.inverse.filelist)"
mdp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"
index=$(csg_get_property cg.inverse.gromacs.grompp.index "index.ndx")
[ -f "$index" ] || die "${0##*/}: grompp index file '$index' not found"
ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
traj="traj.${ext}"

# Prepare and submit simulations
echo "#dist.xvg grofile delta" > dist_comp.d
for i in conf_start*.gro; do
  number=${i#conf_start}
  number=${number%.gro}
  [ -z "$number" ] && die "${0##*/}: Could not fetch number"
  echo Simulation $number
  dir="$(printf sim_%03i $number)"
  mkdir $dir
  cp_from_main_dir $filelist
  mv $i ./$dir/conf.gro
  dist=2
  cp grompp.mdp > $dir/grompp.mdp
  cd $dir
  cp_from_main_dir $filelist
  grompp -n ${index} ${mdp_opts}
  echo -e "${pullgroup0}\n${pullgroup1}" | g_dist -f conf.gro -s topol.tpr -n ${index} -o dist.xvg
  dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
  [ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
  msg "Doing $dir with dist $dist"
  grompp -n ${index} ${mdp_opts}
  do_external run gromacs_pmf
  sleep 5
  cd ..
done

[ -f "$dir/$confout" ] || die "${0##*/}: Gromacs end coordinate '$confout' not found after running mdrun"

cat dist_comp.d | sort -n > dist_comp.d
awk '{if ($4>0.001){print "Oho in step",$1}}' dist_comp.d
