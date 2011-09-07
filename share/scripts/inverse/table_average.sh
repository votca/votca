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
This script creates averages tables and also calculates the error. 

Usage: ${0##*/} [options] table1 table2 table3 ....

Allowed options:
-h, --help                    show this help
-o, --output NANE             output file name
    --cols NUM                Number of columns per file
                              Default: $cols
    --col-y NUM               y-data column
                              Default: $coly
    --col-x NUM               x-data column
                              Default: $colx
    --clean		      Clean intermediate files

Examples:
* ${0##*/} --output CG-CG.dist.new CG-CG*.dist.new 
EOF
}

#default
out=""
cols=3
colx=1
coly=2
clean="no"

shopt -s extglob
while [[ ${1#-} != $1 ]]; do
 if [[ ${1#--} = $1 && -n ${1:2} ]]; then
    #short opt with arguments here: o
    if [[ ${1#-[o]} != ${1} ]]; then
       set -- "${1:0:2}" "${1:2}" "${@:2}"
    else
       set -- "${1:0:2}" "-${1:2}" "${@:2}"
    fi
 fi
 case $1 in
   -o | --output)
    out="$2"
    shift 2;;
   --cols)
    is_int "$2" || die "${0##*/}: argument after --cols should be an integer"
    cols="$2"
    shift 2;;
   --col-x)
    is_int "$2" || die "${0##*/}: argument after --col-x should be an integer"
    colx="$2"
    shift 2;;
   --col-y)
    is_int "$2" || die "${0##*/}: argument after --col-y should be an integer"
    coly="$2"
    shift 2;;
   --clean)
    clean="yes"
    shift;;
   -h | --help)
    show_help
    exit 0;;
  *)
   die "Unknown option '$1'";;
 esac
done
### end parsing options

[[ -z $1 || -z $2 ]] && die "${0##*/}: Missing arguments"
[[ -z "$out" ]] && die "${0##*/}: No output file given add --output option"

l=""
c=0
tables=()
for f in "$@"; do
  [[ -f $f ]] || die "Could not find table $f"
  t=$(critical mktemp "$f.XXXX")
  critical sed '/^[#@]/d' $f > $t
  [[ -z $l ]] && l=$(critical sed -n '$=' $t)
  [[ $l -eq $(critical sed -n '$=' $t) ]] || die "Number of lines (after comments have been striped) mismatches in $f from $1"
  tables[$c]="$t"
  ((c++))
done
t=$(critical mktemp "table_all.XXXX")
paste "${tables[@]}" > "${t}"
awk -v c1="$colx" -v c2="$coly" -v s="$cols" '
func isnum(x){return(x==x+0)}
{
  sum=0;
  sum2=0;
  c=0;
  for (i=0;i<NF;i+=s) {
    sum+=$(i+c2);
    sum2+=$(i+c2)*$(i+c2);
    c++;
    if ($(c1) != $(c1+i)) {
      print "x-value missmatch",$(c1),"vs. ", $(c1+i), " in line",NR > /dev/stderr;
      exit 1;
    }
  }
  flag="u"
  if (isnum(sum)&&isnum(sum2)) { flag="i" }
  print $1,sum/c,sqrt((sum2-sum*sum/c)/(c*(c-1))),flag;
}' $t > $out || die "${0##*/}: averaging with awk failed"

if [[ $clean = "yes" ]]; then
  rm -f "${tables[@]}" "$t"
fi
