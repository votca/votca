#! /bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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
This script calcs the density for gromacs for the AdResS therm force

Usage: ${0##*/}
EOF
   exit 0
fi

dt=$(get_from_mdp dt "grompp.mdp")
equi_time="$(csg_get_property cg.inverse.gromacs.equi_time 0)"
first_frame="$(csg_get_property cg.inverse.gromacs.first_frame 0)"
nsteps=$(get_from_mdp nsteps "grompp.mdp")

name=$(csg_get_interaction_property name)


mdp="$(csg_get_property cg.inverse.gromacs.mdp "grompp.mdp")"
[ -f "$mdp" ] || die "${0##*/}: gromacs mdp file '$mdp' not found"
adress_type=$(get_from_mdp adress_type "$mdp")

echo "Adress type: $adress_type"


begin="$(awk -v dt=$dt -v frames=$first_frame -v eqtime=$equi_time 'BEGIN{print (eqtime > dt*frames ? eqtime : dt*frames) }')"

#TODO rename that to ALL or "*"
densigroup="$(csg_get_interaction_property inverse.gromacs.density_group "UNDEF")"
g_densopt="$(csg_get_property --allow-empty cg.inverse.gromacs.g_density_options)"

#TODO check this
if [ $densigroup = "UNDEF" ]; then
  index_sel="$name"
  echo "Calculating density for $name"
else
  index_sel="$densigroup"
fi

if [ $adress_type = "sphere" ]; then
  #TODO hardcoded values
  critical csg_spheredens --trj traj.trr --top topol.tpr --bin 0.01 --out "dens.$name.xvg" --begin "${begin}" --molname "$densigroup"
else
  #TODO hardcoded stuff
  dens_prog="g_density -n index.ndx -d x $g_densopt"
  echo "Running $dens_prog"
  echo -e "$index_sel" | critical $dens_prog -b "${begin}" -o "dens.$name.xvg"
fi
