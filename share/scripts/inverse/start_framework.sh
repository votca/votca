#! /bin/bash
#
# Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
this script starts the iterative framework and performs a lot of tests
EOF
  exit 0
fi

#for now, we will replace this function later
die(){ echo "$*" >&2; exit 1; }

#check for VOTCASHARE
[[ -n ${VOTCASHARE} ]] || die "Error: VOTCASHARE not definded"
[[ -d ${VOTCASHARE} ]] || die "VOTCASHARE '$VOTCASHARE' is not a dir"
[[ -d ${VOTCASHARE}/scripts/inverse ]] || die "\$VOTCASHARE/scripts/inverse is not found. Is VOTCASHARE set corectly?"
[[ -f ${VOTCASHARE}/scripts/inverse/inverse.sh ]] || die "Could not find inverse.sh, \$VOTCASHARE/scripts/inverse seem to point to the wrong place!"
[[ -f ${VOTCASHARE}/scripts/inverse/functions_common.sh ]] || die "Could not find default common framework functions (functions_common.sh)"
source "${VOTCASHARE}/scripts/inverse/functions_common.sh" || exit 1 

#this is needed by die later
export CSG_MASTER_PID="$$"

export CSG_MAINDIR="$PWD"

if [[ -n ${VOTCA_CSG_DEFAULTS} ]]; then
  [[ -f ${VOTCA_CSG_DEFAULTS} ]] || die "Could not find ${VOTCA_CSG_DEFAULTS}! Is VOTCA_CSG_DEFAULTS set corectly?"
else
  export VOTCA_CSG_DEFAULTS="${VOTCASHARE}/xml/csg_defaults.xml"
  [[ -f ${VOTCA_CSG_DEFAULTS} ]] || die "Could not find ${VOTCA_CSG_DEFAULTS}! Is VOTCASHARE ($VOTCASHARE) set corectly? Hint: When overwriting VOTCASHARE you need to overwrite VOTCA_CSG_DEFAULTS as well."
fi

#do no overwrite CSGSHARE stuff set by user from the outside
add_to_csgshare --at-the-end "${VOTCASHARE}/scripts/inverse"
