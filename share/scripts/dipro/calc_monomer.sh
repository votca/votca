#!/bin/bash
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

###############################################################################
# script to configure and run TURBOMOLE or GAUSSIAN for monomers
###############################################################################

#========================================================================
#  adjusting amount of memory -> RI approximation
#========================================================================
#max. amount of memory currenty 750 MB
adapt_mem(){
  memmax=750
  awk -v mem=$memmax  '
   /ricore/ {printf("$ricore %d\n",mem);next}
   /end/ {print "$end";exit}
   {print $0}
  ' control >temp
  mv temp control
}
#========================================================================

#========================================================================
# running define to prepare turbomole input-files; coord has to be there
#========================================================================
run_define(){

define <<!


a coord
*
no
b all $1
*
eht
y
0
y
dft
on
func
$2
grid
m3
*
ri
on
m 300
*
scf
conv
7
iter
200

marij

q
!

}
#========================================================================


#========================================================================
#   MAIN script execution
#========================================================================


if [ "$1" == "T" ]; then 

   # check if Turbomole environment is set up, if not: abort
    if [ -z "$TURBODIR" ]; then
	echo "Turbomole environment not set! Aborting!"
	exit
    fi

    #single point monomer 
    input=$2
    if [ -n "$input" ]; then
	func=${input%/*}    
	basis=${input##*/}
    else
	basis="def-TZVP"
	func="pbe"
    fi

    XYZ=`ls *.xyz`
    babel -ixyz $XYZ -otmol coord
    run_define $basis $func

    adapt_mem
    kdg scfdump
    kdg end
    kdg scfdamp
    echo "\$scfdamp   start=0.800  step=0.01  min=0.01" >> control
    echo "\$end" >> control
    ridft > ridft.out
    
elif [ "$1" == "G" ]; then

    if [ -z "$g09root" ]; then
	echo "Gaussian 09 environment not set! Aborting!"
	exit
    fi

    XYZ=`ls *.xyz` 
    if [ -n "$2" ]; then
	input=$2
    else
	input="pbepbe/6-311G**"
    fi

    echo "%nproc=1" > monomer.com
    echo "# pop=minimal $input nosymm punch=mo" >> monomer.com
    echo "" >> monomer.com
    echo "monomer generated from calc_monomer" >> monomer.com
    echo "" >> monomer.com
    echo "0 1" >> monomer.com
    sed -e "1,2d" $XYZ >> monomer.com
    echo "" >> monomer.com
    mkdir -p temp
    g09 monomer
    rm -rf temp

fi
























