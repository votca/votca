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

show_help () {
  cat <<EOF
${0##*/}, version %version%
This script creates a table with grid min:step:max for the a functional form

Usage: ${0##*/} [options] outfile

Allowed options:
-h, --help                    show this help
-o, --output NANE             Output file name
    --grid  XX:XX:XX          Output grid of the table
    --var X=Y                 Set a variable used in the function
    --fct FCT                 functional form of the table
    --header XXX              Extra headerlines
    --gnuplot CMD             Gnuplot command to use
                              Default: $gnuplot
    --clean		      Clean intermediate files

Used external packages: gnuplot

Examples:
* ${0##*/} --output CG-CG.dist.new CG-CG*.dist.new 
EOF
}

#default
output=""
vars=()
values=()
grid=()
fct=""
gnuplot=gnuplot
clean="no"
header="$(get_table_comment | sed -e 's/^/#/')\n#\n#Plot script:"

# parse arguments
shopt -s extglob
while [[ ${1} = -* ]]; do
  if [[ ${1} = --*=* ]]; then # case --xx=yy
    set -- "${1%%=*}" "${1#*=}" "${@:2}" # --xx=yy to --xx yy
  elif [[ ${1} = -[^-]?* ]]; then # case -xy split
    if [[ ${1} = -[o]* ]]; then #short opts with arguments
       set -- "${1:0:2}" "${1:2}" "${@:2}" # -xy to -x y
    else #short opts without arguments
       set -- "${1:0:2}" "-${1:2}" "${@:2}" # -xy to -x -y
    fi
 fi
 case $1 in
   -h | --help)
    show_help
    exit 0;;
   -o | --output)
    output="$2"
    shift 2;;
   --fct)
    fct="$2"
    shift 2;;
   --header)
    header+="$header\n$2"
    shift 2;;
   --grid)
    grid[0]="$2"
    [[ $grid =~ ^(.*):(.*):(.*)$ ]] || die "${0##*/}: Agrument after --grid should have the form XX:XX:XX"
    for i in 1 2 3; do
      [[ -z ${BASH_REMATCH[$i]} ]] && die "${0##*/}: part nr $i of XX:XX:XX was empty"
      is_num "${BASH_REMATCH[$i]}" || die "${0##*/}: part nr $i of XX:XX:XX should be a number"
      grid[$i]="${BASH_REMATCH[$i]}"
    done
    shift 2;;
   --clean)
    clean="yes"
    shift;;
   --var)
    [[ ${2//[^=]} = "=" ]] || die "${0##*/}: Agrument after --var should have the form XX=YY"
    [[ -n $vars ]] && vars+="$vars\n$2" || vars="$2"
    shift 2;;
   --debug)
    set -x
    shift ;;
  *)
   die "Unknown option '$1'"
   exit 1;;
 esac
done

for i in grid output fct; do
 [[ -z ${!i} ]] && die "${0##*/}: Missing arguments --$i"
done

len=$(csg_calc "${grid[3]}" - "${grid[1]}")
samples="$(csg_calc "$len" / "${grid[2]}")"
samples="$(to_int "$samples")"
((samples++))

[ -n "$(type -p $gnuplot)" ] || die "${0##*/}: gnuplot binary '$gnuplot' not found"

tmpfile="$(critical mktemp table.gp.XXX)"
tmpfile2="$(critical mktemp table.plot.XXX)"
echo -e "$header" > "$tmpfile"
[[ -n $vars ]] && echo -e "$vars" >> "$tmpfile"
echo "set samples $samples" >> "$tmpfile"
echo "set table '$tmpfile2'" >> "$tmpfile"
echo "plot [${grid[1]}:${grid[3]}] $fct" >> "$tmpfile"
critical $gnuplot "$tmpfile"

critical sed -e 's/^#*/#/' "$tmpfile" > "$output"
echo -e "#\n# Gnuplot output:" >> "$output"
critical sed -e '/^[[:space:]]*$/d' "$tmpfile2" >> "$output"

if [[ $clean = "yes" ]]; then
  rm -f "$tmpfile" "$tmpfile2"
fi
