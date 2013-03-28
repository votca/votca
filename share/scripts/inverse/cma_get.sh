#! /bin/bash
#
# Copyright 2009-2012 The VOTCA Development Team (http://www.votca.org)
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

if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
Will download cma python module and put it in CSGSHARE

Usage: ${0##*/}
EOF
   exit 0
fi

[[ -n "$(type -p wget)" ]] || die "${0##*/}: Could not find wget"
msg "We will now go to http://www.lri.fr/~hansen/cmaes_inmatlab.html and get the cma python module"
msg "Please consider sending Nikolaus Hansen an email telling him that you are using his code"
msg "within VOTCA. hansen@lri.fr, it will only take seconds!"
sleep 10
critical wget -O "${CSGSHARE}/cma.py" "http://www.lri.fr/~hansen/cma.py"
critical chmod 755 "${CSGSHARE}/cma.py"
msg "Done"
