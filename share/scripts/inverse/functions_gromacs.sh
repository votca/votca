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

if [ "$1" = "--help" ]; then
  #we add \$GMXDATA in USES, because gromacs will need it
cat <<EOF
${0##*/}, version %version%
Useful functions for gromacs 4.0

NEEDS:

USES: sed die \$GMXDATA check_deps

PROVIDES: get_from_mdp cufoff_check
EOF
  exit 0
fi

check_deps $0

get_from_mdp() {
  local res
  [[ -n "$2" ]] || die "get_from_mdp: Missing argument (what file)"
  [[ -f "$2" ]] || die "get_from_mdp: Could not read file '$2'"
  #1. strip comments
  #2. get important line
  #3. remove leading and tailing spaces
  res="$(sed -e '/^[[:space:]]*;/d' -e 's#;.*$##' "$2" | \
        sed -n -e "s#^[[:space:]]*$1[[:space:]]*=[[:space:]]*\(.*\)[[:space:]]*\$#\1#p" | \
	sed -e 's#^[[:space:]]*##' -e 's#[[:space:]]*$##')" || \
    die "get_from_mdp: sed failed"
  [[ -z "$res" ]] && [ -z "$3" ] && die "get_from_mdp: could not fetch $1 from $2"
  [ -n "$res" ] && echo "$res" || echo "$3"
}
export -f get_from_mdp

check_cutoff() {
  local max rvdw res cutoff_check
  [[ -n "$1" ]] || die "check_cutoff: Missing argument (interaction name)"
  cutoff_check=$(csg_get_property cg.inverse.gromacs.cutoff_check "yes")
  [ ${cutoff_check} = "no" ] && return 0
  max="$(csg_get_interaction_property max)"
  rvdw="$(get_from_mdp rvdw "$1" 1)"
  res="$(awk -v max="$max" -v rvdw="$rvdw" 'BEGIN{ print (max>rvdw)?1:0 }')" || die "check_cutoff: awk failed"
  [ "$res" != "0" ] && die "Error in interaction '$bondname': rvdw ($rvdw) in $1 is smaller than max ($max)\n\
To ignore this check set cg.inverse.gromacs.cutoff_check to 'no'"
  return "$res"
}
export -f check_cutoff
