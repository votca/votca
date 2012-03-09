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


if [[ $1 = "--help" ]]; then
cat <<EOF
${0##*/}, version %version%
Checks whether csg_reupdate is converged by looking for file 'converged', if it exists then 
creates a file 'stop' which is seen by inverse script and program is stopped.

Usage: ${0##*/}
EOF
   exit 0
fi

# TO DO

#out='converged'
#stp='stop'
#if [[ -f $out ]]; then
#   mv $out $stp
#fi

#if [[ -f notsympos ]]; then
#   msg --color red "Hessian not a positive definite"
#   if [ "$(csg_get_property cg.inverse.re.dosteepest)" = "no" ]; then
#    msg --color red "User decided not to take steepest descent"
#    touch $stp
#   else
#    msg --color red "User is ok with taking steepest descent"
#   fi  
#fi
