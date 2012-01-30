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
This script implements the preparation of the relative entropy method iteration

Usage: ${0##*/}
EOF
   exit 0
fi


#get initial parameters from main dir and make it current parameters
for_all non-bonded 'cp_from_main_dir --rename $(csg_get_interaction_property name).param.init $(csg_get_interaction_property name).param.cur'

# copy all mapping xml file from the main directory
# need to do this more elegantly
cp_from_main_dir map*.xml 

# copy coarse-grained initial configurations
cp_from_main_dir --rename conf.gro confout.gro

# generate reference AA ensemble CG-CG histograms and 
# also generate potential tables from the initial parameter guess

msg --color green "Generating AA ref histograms, initial CG topology, and potential tables from the initial guess"

tpr="$(get_main_dir)/$(csg_get_property cg.ref.top)"
[[ -f $tpr ]] || die "${0##*/}: please provide ref AA topology file '$tpr' "

trr="$(get_main_dir)/$(csg_get_property cg.ref.trj)"
[[ -f $trr ]] || die "${0##*/}: please provide ref AA trajectory file '$trr' "

ext=$(csg_get_property cg.inverse.gromacs.traj_type)
[[ $ext = trr ]] || die "Method re needs a trr trajetory (to read forces)"

mapping="$(csg_get_property cg.inverse.map)"
opts="$(csg_get_property cg.inverse.re.csg_reupdate.opts)"

# for now I am providing initial step cg conf
# this wont be needed if i could find better way to do this
# working on it

critical csg_reupdate --top ${tpr} --trj ${trr} --options $CSGXMLFILE --genref true --cg $(csg_get_property cg.map) ${opts}
