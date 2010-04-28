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

die "3 parameters are nessary\n" if ($#ARGV<2);

use CsgFunctions;
use SimplexFunctions;

my $outfile="$ARGV[0]";
my $simplex_table="$ARGV[1]";
my $p_line_nr="$ARGV[2]";

my $min=csg_get_property("cg.non-bonded.min");
my $max=csg_get_property("cg.non-bonded.max");
my $step=csg_get_property("cg.non-bonded.step");

my @r;
my @dummy;
my @flag;

open(TMP, ">tmp");
print TMP "$min 0\n$max 0";

my @args=("bash","-c","csg_resample --in tmp --out grid --grid $min:$step:$max");
system(@args);

(readin_table("grid",@r,@dummy,@flag)) || die "$progname: error at readin_table\n";

close(TMP);

my @ftar;
my @sig;
my @eps;
my @flag_simplex;

(readin_simplex_table($simplex_table,@ftar,@sig,@eps,@flag_simplex)) || die "$progname: error at readin_simplex_table\n";

my @pot;
for (my $i=0;$i<=$#r;$i++){
    # Avoid undefined potential at r=0
    if ($r[$i]>1e-10) {
        $pot[$i]=calc_func("$r[$i]","$sig[$p_line_nr]","$eps[$p_line_nr]");
        $flag[$i]="i";
    }
    else {
      $pot[$i]="0";
      $flag[$i]="u";
    }
    # Avoid gmx segmentation fault for large pot
    if ($pot[$i]>=1e6) {
        $pot[$i]=1e6;
    }
}

# Find index at the cutoff
my $i_cut=$#r;

# Shift potential so that it is zero at the cutoff
for (my $i=0;$i<=$i_cut;$i++){
   $pot[$i]-=$pot[$i_cut];
}

$flag_simplex[$p_line_nr]="active";

saveto_table($outfile,@r,@pot,@flag) || die "$progname: error at saveto_table\n";
saveto_simplex_table("simplex.new",@ftar,@sig,@eps,@flag_simplex) || die "$progname: error at saveto_simplex_table\n";