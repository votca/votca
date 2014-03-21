#! /bin/bash
#
# Copyright 2009-2014 The VOTCA Development Team (http://www.votca.org)
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
This script converts all potentials to the format needed by dlpoly

Usage: ${0##*/}
EOF
   exit 0
fi

for i in TABLE TABBND TABANG TABDIH; do
  if [[ -f "${i}" ]]; then
    echo "Storing ${i} ..."
    if [[ -f "${i}.prv" ]]; then
      echo "Overwriting ${i}.prv ..."
    fi
    critical mv -vf "${i}" "${i}.prv"
  fi
done

for_all "non-bonded bonded" do_external convert_potential dlpoly '$(csg_get_interaction_property name).pot.cur' '$(csg_get_interaction_property name).pot.dlpoly'

for i in TAB*prv; do
  if [[ -s "${i}" ]]; then
    fname="$(basename ${i} .prv)"
    if [[ -s "${fname}" ]]; then
      echo "New ${fname} has been created based on the latest distribution"
    else
      echo "Restoring ${fname} ..."
      critical mv -vf "${i}" "${fname}"
    fi
  fi
done
