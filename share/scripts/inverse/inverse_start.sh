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
this script just contain all startup check for inverse.sh to
and makes inverse.sh more readable

PROVIDES: \$CSGXMLFILE \$CSGINVERSE \$SOURCE_WRAPPER \$CSGSCRIPTDIR \$CSGLOG \$CSGRESTART

USES: \$CSGSHARE csg_property die csg_get_property rm int_check

NEEDS: cg.inverse.scriptdir cg.inverse.log_file cg.inverse.restart_file
EOF
  exit 0
fi

#no check deps $0 here
#because this bootstrap everything

#for now, we will replace this function later
die(){ echo "$*" >&2; exit 1; }

[[ -n "$1" ]] || die "Error: Missing xml file"

if [ -f "./$1" ]; then
  export CSGXMLFILE="${PWD}/${1}"
else
  die "Error: file '$1' could not read, needed for \$CSGXMLFILE"
fi

#check for CSGSHARE
[[ -n "$CSGSHARE" ]] || die "Error: CSGSHARE not definded"
[[ -d "$CSGSHARE" ]] || die "CSGSHARE '$CSGSHARE' is not a dir"

CSGINVERSE="${CSGSHARE}/scripts/inverse"
[[ -d "$CSGINVERSE" ]] || die "CSGSHARE/scripts/inverse is not a dir"
export CSGINVERSE
export PERL5LIB="$CSGINVERSE:$PERL5LIB"

#we need csg_property
[[ -n "$(type -p csg_property)" ]] || die "Error: csg_property not found, check your PATH"

#find source_wrapper.pl
SOURCE_WRAPPER="${CSGINVERSE}/source_wrapper.pl"
[[ -x "${SOURCE_WRAPPER}" ]] || die "Could not find source_wrapper.pl"
export SOURCE_WRAPPER

function_file=$($SOURCE_WRAPPER functions common) || die "$SOURCE_WRAPPER functions common failed"
#die() is overwritten here
source ${function_file} || exit 1
unset function_file

CSGSCRIPTDIR="$(csg_get_property --allow-empty cg.inverse.scriptdir)"
#scriptdir maybe contains $PWD or something
if [ -n "$CSGSCRIPTDIR" ]; then
  eval CSGSCRIPTDIR=$CSGSCRIPTDIR
  CSGSCRIPTDIR="$(cd $CSGSCRIPTDIR;pwd)"
  [[ -d "$CSGSCRIPTDIR" ]] || die "CSGSCRIPTDIR '$CSGSCRIPTDIR' is not a dir"
  export CSGSCRIPTDIR
  export PERL5LIB="$CSGSCRIPTDIR:$PERL5LIB"
fi

CSGLOG="$(csg_get_property cg.inverse.log_file)"
CSGLOG="$PWD/${CSGLOG##*/}"
export CSGLOG

#define $CSGRESTART
CSGRESTART="$(csg_get_property cg.inverse.restart_file)"
CSGRESTART="${CSGRESTART##*/}"
export CSGRESTART

export CSG_MAINDIR="$PWD"

#stuff for options
[ -z "$do_iterations" ] || int_check "$do_iterations" "inverse.sh: --do-iterations need a number as agrument"
[ -z "$wall_time" ] || int_check "$wall_time" "inverse.sh: --wall-time need a number as agrument"
if [ "$clean" = "yes" ]; then
  echo -e "So, you want to clean?\n"
  echo "We will remove:"
  files="$(ls -d done ${CSGRESTART} ${CSGLOG##$PWD/} step_* *~ 2>/dev/null)"
  if [ -z "$files" ]; then
    echo "Nothing to clean"
  else
    echo $files
    echo -e "\nCTRL-C to stop it"
    for ((i=10;i>0;i--)); do
      echo -n "$i "
      sleep 1
    done
    rm -rf $files
    echo -e "\n\nDone, hope you are happy now"
  fi
  exit 0
fi
