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
Changes a simplex state according to the current state and current table of parameter using the Nelder–Mead method or downhill simplex.

Usage: ${0##*/} currrentstate currenttable newstate
EOF
   exit 0
fi

#move this to separate file
print_state() {
  [[ -z $1 ]] && die "${FUNCNAME[0]}: missing argument"
  [[ $1 = "Contraction" && -z $2 ]] && die "${FUNCNAME[0]}: Contraction needs an extra arguments"
  echo "State: $1"
  [[ $1 = "Contraction" ]] && echo "Saved: $2"
}

get_state() {
  local state
  [[ -z $1 ]] && die "${FUNCNAME[0]}: missing argument"
  [[ -f $1 ]] || die "${FUNCNAME[0]}: Could not read '$1'"
  state=$(critical awk '/^State:/{print $2}' "$1")
  [[ -z "$state" ]] && die "${FUNCNAME[0]}: Could not get state from state file '$1'"
  extra=$(critical awk '/^Saved:/{print $2}' "$1")
  echo "$state"
  echo "$extra"
}

min_value() {
  local max
  for i in $*; do
    is_num "$i" || die "${FUNCNAME[0]}: $i is not a number"
    [[ -z $min ]] && min="$i"
    csg_calc "$i" "<" "$min" && min="$i"
  done
  echo "$min"
}

max_value() {
  local max
  for i in $*; do
    is_num "$i" || die "${FUNCNAME[0]}: $i is not a number"
    [[ -z $max ]] && max="$i"
    csg_calc "$i" ">" "$max" && max="$i"
  done
  echo "$max"
}

get_highest_value() {
  local values
  [[ -z $1 ]] && die "${FUNCNAME[0]}: missing argument"
  [[ -f $1 ]] || die "${FUNCNAME[0]}: Could read '$1'"
  values=$(critical awk -F '#' '/^[^#@]/{print $NF;}' $1 | critical awk '{print $1}')
  [[ -z $values ]] && die "${FUNCNAME[0]}: Could not get values from table files '$1'"
  max "$values"
}

get_try_value() {
  local value
  [[ -z $1 ]] && die "${FUNCNAME[0]}: missing argument"
  [[ -f $1 ]] || die "${FUNCNAME[0]}: Could read '$1'"
  value=$(critical awk -F '#' '/try$/{print $NF;}' $1 | critical awk '{print $1}')
  is_num "$value" || die "${FUNCNAME[0]}: Could not get try value from table files '$1'"
  echo "$value"
}

[[ -z $1 || -z $2 || -z $3 ]] && die "${0##*/}: Missing argument"
state="$1"
state_new="$3"
table="$2"

[[ -f $state ]] || die "${0##*/}: Could not find current simplex state file '$state'"
[[ -f $table ]] || die "${0##*/}: Cound not find current simplex table file '$state'"

this_state=( $(get_state $state) )
case "${this_state[0]}" in
  None)
    next_state="Reflection";;
  Reflection)
 
}

[[ -z $1 || -z $2 || -z $3 ]] && die "${0##*/}: Missing argument"
state="$1"
state_new="$3"
table="$2"

[[ -f $state ]] || die "${0##*/}: Could not find current simplex state file '$state'"
[[ -f $table ]] || die "${0##*/}: Cound not find current simplex table file '$state'"

this_state=( $(get_state $state) )
try=$(get_try_value $table)
case "${this_state[0]}" in
  None)
    next_state="Reflection";;
  Reflection)
    if csg_calc "$try" > "$(get_highest_value $table)"; then
      next_state="Contraction"
    elif csg_calc "$try" < "$(get_lowest_value $table)"; then
      next_state="Expansion"
    else
      next_state="Reflection"
    fi;;
  Expansion)
  Contraction)
  Reduction)
  *)
  die "${0##*/}: I don't know what to in case of state '$this_state'";;
esac
