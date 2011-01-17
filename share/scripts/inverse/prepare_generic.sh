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
This script initializes potentials in a generic way

Usage: ${0##*/}

USES:  csg_get_property for_all do_external check_deps

NEEDS: cg.inverse.method cg.inverse.program
EOF
   exit 0
fi

check_deps "$0"

sim_prog="$(csg_get_property cg.inverse.program)"
method="$(csg_get_property cg.inverse.method)"

for_all non-bonded do_external prepare_single $method

#cp confout.gro and so on
do_external prepare_generic $sim_prog
