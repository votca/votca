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
This script implemtents the function initialize
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES: do_external csg_get_interaction_property log run_or_exit csg_resample log

NEEDS: name min max step
EOF
  exit 0
fi

check_deps "$0"

name=$(csg_get_interaction_property name)
if [ -f ../${name}.pot.in ]; then
  msg "Using given table ${name}.pot.in for ${name}"
  min=$(csg_get_interaction_property min )
  max=$(csg_get_interaction_property max )
  step=$(csg_get_interaction_property step )
  run_or_exit csg_resample --in ../${name}.pot.in --out ${name}.pot.new --grid ${min}:${step}:${max}
else
  # RDF_to_POT.pl just does log g(r) + extrapolation
  msg "Using intial guess from RDF for ${name}"
  run_or_exit do_external rdf pot ${name}.dist.tgt ${name}.pot.new
fi

