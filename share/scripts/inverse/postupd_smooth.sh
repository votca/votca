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
This script implemtents smoothing of the potential update (.dpot)

Usage: ${0##*/}

USES:  die csg_get_interaction_property mktemp do_external cp log run_or_exit

NEEDS: name inverse.post_update_options.smooth.iterations
EOF
   exit 0
fi

check_deps "$0"

name=$(csg_get_interaction_property name)
tmpfile=$(mktemp ${name}.XXX) || die "mktemp failed"
iterations=$(csg_get_interaction_property inverse.post_update_options.smooth.iterations 1)  

run_or_exit cp ${name}.dpot.cur $tmpfile
log "doing $iterations smoothing iterations"

for((i=0;i<$iterations;i++)); do
  run_or_exit do_external table smooth $tmpfile ${name}.dpot.new
  run_or_exit cp ${name}.dpot.new $tmpfile
done

