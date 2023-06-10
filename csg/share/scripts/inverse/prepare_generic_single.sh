#! /bin/bash
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

show_help () {
  cat <<EOF
${0##*/}, version %version%
This script prepares the potential in step 0, using pot.in or by resampling and
inverting the target distribution

Use --use-table or --use-bi to enforce the method. Otherwise it will use
.pot.in if present and BI if not.

If using --table-overwrite it will overwrite the initial guess with a table
if one is present.

Usage: ${0##*/} [--help] [--use-table|--use-bi]"
EOF
}

use_table=false
use_bi=false
table_overwrite=false
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
  --use-table)
    use_table=true
    shift  # past argument
    ;;
  --use-bi)
    use_bi=true
    shift  # past argument
    ;;
  --table-overwrite)
    table_overwrite=true
    shift  # past argument
    ;;
  --help)
    show_help
    exit 0
    ;;
  *)
    die "unknown argument $key"
    ;;
  esac
done

if [[ (${use_table} == true && ${use_bi} == true)
  || (${use_table} == true && ${table_overwrite} == true)
  || (${use_bi} == true && ${table_overwrite} == true) ]]; then
  die "use either --use-table or --use-bi or --table-overwrite"
fi

name=$(csg_get_interaction_property name)
min=$(csg_get_interaction_property min )
max=$(csg_get_interaction_property max )
step=$(csg_get_interaction_property step )
comment="$(get_table_comment)"
main_dir=$(get_main_dir)
bondtype="$(csg_get_interaction_property bondtype)"
output="${name}.pot.new"
table_overwrite="$(csg_get_property cg.inverse.initial_guess.table_overwrite)"

table_present=false
if [[ -f ${main_dir}/${name}.pot.in ]]; then
  table_present=true
fi

# if --use-bi was used
if [[ $use_bi == true ]]; then
  if [[ $table_present == true && $table_overwrite == true ]]; then
    do_external initial_guess_single table
  elif [[ $table_present == true && $table_overwrite == false ]]; then
    msg --color blue "###########################################################################################"
    msg --color blue "# WARNING there is a table ${name}.pot.in present, but cg.inverse.initial_guess.method=bi #"
    msg --color blue "# and cg.inverse.initial_guess.table_overwrite=false. Using BI as initial guess           #"
    msg --color blue "###########################################################################################"
    do_external initial_guess_single bi
  else
    do_external initial_guess_single bi
  fi
  # if --use-table was used
elif [[ $use_table == true ]]; then
  do_external initial_guess_single table
  # if --table-overwrite was used
elif [[ $table_overwrite == true ]]; then
  # use table if table exists
  if [[ $table_present == true ]]; then
    msg "there is a table ${name}.pot.in present, gonna use it"
    do_external initial_guess_single table
  fi
  # if no table present will do nothing
else
  # this is the old default behaviour
  # use table if table exists
  if [[ $table_present == true ]]; then
    msg "there is a table ${name}.pot.in present, gonna use it"
    do_external initial_guess_single table
  # use Boltzmann inversion otherwise
  else
    do_external initial_guess_single bi
  fi
fi
