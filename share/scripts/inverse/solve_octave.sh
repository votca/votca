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
This solves linear equation system from imc using octave

Usage: ${0##*/} <group> <outfile>

USES:  die sed octave rm run_or_exit cat_external

NEEDS:
EOF
   exit 0
fi

check_deps "$0"

[[ -n "$2" ]] || die "${0##*/}: Missing arguments"

# initialize & run the octave file
cat_external solve octave | sed -e "s/\$name_out/$2/"  -e "s/\$name/$1/" > solve_$1.octave || die "${0##*/}: sed failed"
run_or_exit octave solve_$1.octave
#rm -f solve_$1.octave

[[ -f "$2" ]] || die "Octave failed"
# temporary compatibility issue
run_or_exit sed -ie 's/NaN/0.0/' $2
run_or_exit sed -ie 's/Inf/0.0/' $2

