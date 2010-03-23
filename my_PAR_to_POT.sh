#! /usr/bin/perl -w
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

use strict;

( my $progname = $0 ) =~ s#^.*/##;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
$progname, %version%
This script calculates the potential from parameters
given in file simplex.new
In addition it does some magic tricks:
- does not crash when calc r(0)
- does not allow values of pot>1e10
- shifts the potential, so that it is zero at the cutoff

Usage: $progname infile outfile

USES: readin_table readin_simplex_table calc_func saveto_table

NEEDS: -

EOF
  exit 0;
}

die "4 parameters are nessary\n" if ($#ARGV<3);

use CsgFunctions;
use SimplexFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
my $simplex_table="$ARGV[2]";
my $a_line_nr="$ARGV[3]";

my @r;
my @rdf;
my @flag;
(readin_table($infile,@r,@rdf,@flag)) || die "$progname: error at readin_table\n";

my @f_target;
my @sig;
my @eps;
(readin_simplex_table("$simplex_table",@f_target,@sig,@eps)) || die "$progname: error at readin_simplex_table\n";

my @pot;
for (my $i=0;$i<=$#r;$i++){
    # Avoid undefined potential at r=0
    if ($r[$i]>1e-10) {
        $pot[$i]=calc_func("$r[$i]","$sig[$a_line_nr]","$eps[$a_line_nr]");
        $flag[$i]="i";
    }
    else {
      $pot[$i]="0";
      $flag[$i]="u";
    }
    # Avoid gmx segmentation fault for large pot
    if ($pot[$i]>=1e10) {
        $pot[$i]=1e10;
    }
}

# Find index at the cutoff
my $i_cut=$#r;

# Shift potential so that it is zero at the cutoff
for (my $i=0;$i<=$i_cut;$i++){
   $pot[$i]-=$pot[$i_cut];
}

saveto_table($outfile,@r,@pot,@flag) || die "$progname: error at save table\n";