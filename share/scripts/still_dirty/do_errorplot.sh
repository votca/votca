#!/bin/bash
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

do_errors() {
  csg_imc --top topol.tpr --trj traj.xtc --cg ../water_cg.xml --options ../cg.xml \
    --do-imc --begin 100 \
    --write-every $1 --do-blocks --nframes $((16*$1))
  ~/src/csg/share/scripts/inverse/errorbars.sh ibm
  cp CG-CG.dpot.err ibm.err.$i
  ~/src/csg/share/scripts/inverse/errorbars.sh imc
  cp CG-CG.dpot.err imc.err.$i
}

for i in $(seq 1 250); do
  do_errors $i
done
