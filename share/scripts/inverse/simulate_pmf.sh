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

pullgroup0=$(csg_get_interaction_property .pmf.pullgroup0)
pullgroup1=$(csg_get_interaction_property .pmf.pullgroup1)
conf_start="start"
min=$(csg_get_interaction_property .pmf.from)
max=$(csg_get_interaction_property .pmf.to)
steps=$(csg_get_interaction_property .pmf.steps)
rate=$(csg_get_interaction_property .pmf.rate)
out=$(csg_get_interaction_property .pmf.out)
sim_time=$(csg_get_interaction_property .pmf.sim_time)
step_time=$(csg_get_interaction_property .pmf.step_time)
filelist="$(csg_get_property --allow-empty cg.inverse.filelist)"
mdp_init="start_in.mdp"
mdp_opts="$(csg_get_property --allow-empty cg.inverse.gromacs.grompp.opts)"
ext=$(csg_get_property cg.inverse.gromacs.traj_type "xtc")
traj="traj.${ext}"

# Generate start_in.mdp
cp_from_main_dir grompp.mdp.template
cat grompp.mdp.template | sed 's/^pull.*$//' | uniq > tmp
sed -e "s/@TIMESTEP@/$dt/" \
    -e "s/@EXCL@//" \
    -e "s/@STEPS@/$steps/" tmp > ${mdp_init}
rm tmp

dt=$(get_from_mdp dt "$mdp_init")

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
      -e "s/@OUT@/$out/" \
      -e "s/@EXCL@//" \
      -e "s/@STEPS@/$steps/" grompp.mdp.template > $dir/grompp.mdp
  cd $dir
  cp_from_main_dir $filelist
  grompp -n index.ndx ${mdp_opts}
  echo -e "${pullgroup0}\n${pullgroup1}" | g_dist -f conf.gro -s topol.tpr -n index.ndx -o dist.xvg
  dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
  [ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
  msg "Doing $dir with dist $dist and time $step_time"
  critical sed -e "s/@DIST@/$dist/" \
      -e "s/@RATE@/0/" \
      -e "s/@TIMESTEP@/$dt/" \
      -e "s/@OUT@/$out/" \
      -e "s/@EXCL@//" \
      -e "s/@STEPS@/$steps/" ../grompp.mdp.template > grompp.mdp
  grompp -n index.ndx ${mdp_opts}
  do_external run gromacs_pmf
  sleep 5
  cd ..
done

# Wait for jobs to finish and resubmit
sleep 10

times=$(seq -w $step_time $step_time $sim_time)
msg "---------------"
for t in $times; do  
  for dir in sim_*; do
    cd $dir
    dist=$(sed '/^[#@]/d' dist.xvg | awk '{print $2}')
    [ -z "$dist" ] && die "${0##*/}: Could not fetch dist"
    background=$(csg_get_property cg.inverse.simulation.background "no")
    sleep_time=$(csg_get_property cg.inverse.simulation.sleep_time "60")
    if [ "$background" == "yes" ]; then
      while [ -z $(find . -maxdepth 1 -name confout.part*) ]; do
        sleep $sleep_time
      done
      if [ "$t" == "$sim_time" ]; then
        msg "$dir is done"
        confout=$(basename $(find . -maxdepth 1 -name confout.part*) )
        traj=$(basename $(find . -maxdepth 1 -name traj.part*) )
        mv topol.tpr topol.pg.${t}.tpr
        mv grompp.mdp grompp.pg.${t}.mdp
        touch sim_done

        # Rerun finished trajectory
        msg "Doing $dir with mdrun -rerun"
        cp_from_main_dir grompp.mdp.template
        critical sed -e "s/@DIST@/$dist/" \
        -e "s/@RATE@/0/" \
        -e "s/@TIMESTEP@/$dt/" \
        -e "s/@OUT@/0/" \
        -e "s/@STEPS@/$steps/" \
        -e "s/@EXCL@/pullgroup0 environment pullgroup1 environment pullgroup0 pullgroup0 pullgroup1 pullgroup1 environment environment/" grompp.mdp.template > grompp.mdp
        rm grompp.mdp.template
        grompp -n index.ndx ${mdp_opts} -o topol.tpr
        do_external run gromacs_pmf
        while [ ! -f "rerun_done" ]; do
          sleep $sleep_time
        done
        rm topol.tpr
        rm sim_done
        rm rerun_done
        rm topol.pg.${t}.tpr
        rm ${traj}
        
        msg "Concatenating files..."
        # Concatenate energies
        eneconv -f ener.part* -o ener.edr
        while [ ! -f "ener.edr" ]; do
          sleep $sleep_time
        done
        rm ener.part*
        # Concatenate energies 2
        eneconv -f ener2.part* -o ener2.edr
        while [ ! -f "ener2.edr" ]; do
          sleep $sleep_time
        done
        rm ener2.part*
        # Concatenate pullf.xvg files
        head -12 pullf.part0001.xvg > pullf.xvg
        cat pullf.part000* | grep -v ^# | grep -v ^@ | uniq >> pullf.xvg
        while [ ! -f "pullf.xvg" ]; do
          sleep $sleep_time
        done
        rm pullf.part*
        # Concatenate pullf2.xvg files
        head -12 pullf2.part0002.xvg > pullf2.xvg
        cat pullf2.part000* | grep -v ^# | grep -v ^@ | uniq >> pullf2.xvg
        while [ ! -f "pullf2.xvg" ]; do
          sleep $sleep_time
        done
        rm pullf2.part*
        cd ..
      else
        confout=$(basename $(find . -maxdepth 1 -name confout.part*) )
        traj=$(basename $(find . -maxdepth 1 -name traj.part*) )
        mv topol.tpr topol.pg.${t}.tpr
        mv grompp.mdp grompp.pg.${t}.mdp
        touch sim_done
        
        # Rerun finished trajectory
        msg "Doing $dir with mdrun -rerun"
        cp_from_main_dir grompp.mdp.template
        critical sed -e "s/@DIST@/$dist/" \
        -e "s/@RATE@/0/" \
        -e "s/@TIMESTEP@/$dt/" \
        -e "s/@OUT@/0/" \
        -e "s/@STEPS@/$steps/" \
        -e "s/@EXCL@/pullgroup0 environment pullgroup1 environment pullgroup0 pullgroup0 pullgroup1 pullgroup1 environment environment/" grompp.mdp.template > grompp.mdp
        rm grompp.mdp.template
        grompp -n index.ndx ${mdp_opts} -o topol.tpr
        do_external run gromacs_pmf
        while [ ! -f "rerun_done" ]; do
          sleep $sleep_time
        done
        rm topol.tpr
        rm sim_done
        rm rerun_done
        
        # Prepare new mdrun
        msg "Doing $dir with dist $dist and time $(($t+$step_time))"
        echo "System" | tpbconv -s topol.pg.${t}.tpr -f ${traj} -e ener.edr -n index.ndx -until $(($t+$step_time)) -o topol.tpr
        while [ ! -f "topol.tpr" ]; do
          sleep $sleep_time
        done
        rm ${traj}
        rm topol.pg.${t}.tpr
        rm conf.gro
        mv $confout conf.gro
        do_external run gromacs_pmf
        cd ..
      fi
    else
      [ -f "$dir/$confout" ] || die "${0##*/}: Gromacs end coordinate '$confout' not found after running mdrun"
    fi
  done
  msg "---------------"
done

cat dist_comp.d | sort -n > dist_comp.d
awk '{if ($4>0.001){print "Oho in step",$1}}' dist_comp.d
