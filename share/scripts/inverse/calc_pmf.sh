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
This script implemtents the function initialize for the PMF calculator

Usage: ${0##*/}

EOF
  exit 0
fi

if [ -f ${CSGSHARE}/scripts/inverse/functions_pmf.sh ]; then
  source ${CSGSHARE}/scripts/inverse/functions_pmf.sh || die "Could not source functions_pmf.sh"
else
  die "Could not find functions_pmf.sh"
fi

last_step="03_simulations"
pullgroup0=$(csg_get_property cg.non-bonded.name | sed 's/-.*$//')
pullgroup1=$(csg_get_property cg.non-bonded.name | sed 's/^.*-//')
kBT=$(csg_get_property cg.inverse.kBT)
conf_start="start"
steps=$(csg_get_property cg.non-bonded.steps)
f_meas=$(csg_get_property cg.non-bonded.f_meas)
out=$((steps/f_meas))

forcefile="forces_${pullgroup0}_${pullgroup1}.d"
potfile="pmf_${pullgroup0}_${pullgroup1}.d"
rdffile="rdf_${pullgroup0}_${pullgroup1}.d"

echo "#dist <force> error flag" > ${forcefile}
sims="$(find ../${last_step} -type d -name "sim_???*" | sort)"
for sim in ${sims}; do
  echo -n "Doing ${sim}"
  dist="$(get_from_mdp pull_init1 ${sim}/grompp.mdp)"
  [ -f "${sim}/pullf.xvg" ] || die "Could not find file ${sim}/pullf.xvg"
  force="$(${CSGSHARE}/scripts/inverse/avg_bl.awk -v col=2 ${sim}/pullf.xvg | awk '/^[^#]/{print $1,$2}')"
  echo "dist: $dist force: $force"
  echo "$dist $force i" >>  ${forcefile}
done

cp ${sim}/grompp.mdp .

echo "Calculating pmf"
# Make forces negative and the integrate to get pot
do_external table_err linearop --with-errors ${forcefile} tmp -1 0
do_external table integrate --with-errors --with-S --kbT $kBT tmp ${potfile}
rm tmp
