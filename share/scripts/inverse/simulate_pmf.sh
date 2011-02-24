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

pullgroup0="pullgroup0"
pullgroup1="pullgroup1"

conf_start="start"
min=$(csg_get_property cg.non-bonded.min)
max=$(csg_get_property cg.non-bonded.max)
steps=$(csg_get_property cg.non-bonded.steps)
dt=$(csg_get_property cg.non-bonded.dt)
rate=$(csg_get_property cg.non-bonded.rate)
f_meas=$(csg_get_property cg.non-bonded.f_meas)
out=$((steps/f_meas))
filelist="$(csg_get_property --allow-empty cg.inverse.filelist)"

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
  critical sed -e "s/@DIST@/$dist/" \
      -e "s/@RATE@/0/" \
      -e "s/@TIMESTEP@/$dt/" \
      -e "s/@OUT@/0/" \
      -e "s/@PULL_OUT@/$out/" \
      -e "s/@STEPS@/$steps/" grompp.mdp.template > $dir/grompp.mdp
  cd $dir
  cp_from_main_dir $filelist
  grompp -n index.ndx
  echo -e "pullgroup0\npullgroup1" | g_dist -f conf.gro -s topol.tpr -n index.ndx -o dist.xvg
  dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
  [ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
  msg "Dist is $dist"
  critical sed -e "s/@DIST@/$dist/" \
      -e "s/@RATE@/0/" \
      -e "s/@TIMESTEP@/$dt/" \
      -e "s/@OUT@/0/" \
      -e "s/@PULL_OUT@/$out/" \
      -e "s/@STEPS@/$steps/" ../grompp.mdp.template > grompp.mdp

  grompp -n index.ndx
  do_external run gromacs_pmf
  cd ..
done

# Wait for jobs to finish
sleep 10
for dir in sim_*; do
  dir="$(printf sim_%03i $number)"
  confout="$(csg_get_property cg.inverse.gromacs.conf_out "confout.gro")"
  background=$(csg_get_property cg.inverse.parallel.background "no")
  sleep_time=$(csg_get_property cg.inverse.parallel.sleep_time "60")
  if [ "$background" == "yes" ]; then
    while [ ! -f "$dir/$confout" ]; do
      sleep $sleep_time
    done
  else
    ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
    traj="traj.${ext}"
    [ -f "$dir/$confout" ] || die "${0##*/}: Gromacs end coordinate '$confout' not found after running mdrun"
  fi
done

cat dist_comp.d | sort -n > dist_comp.d
awk '{if ($4>0.001){print "Oho in step",$1}}' dist_comp.d
