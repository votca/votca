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
This script calculates the 12-6 Lennard Jones potential
for a given parameter set.

Usage: $progname infile outfile tmp param_N p_line_nr

USES: readin_table csg_get_property csg_resample saveto_table

NEEDS: min max step
EOF
  exit 0;
}

# Give number of parameters
if (defined($ARGV[0])&&("$ARGV[0]" eq "--nparams")){
  print "2\n";
  exit 0;
}

die "5 parameters are nessary\n" if ($#ARGV<4);

use CsgFunctions;
use SimplexFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
my $tmp="$ARGV[2]";
my $param_N="$ARGV[3]";
my $p_line_nr="$ARGV[4]";

my $min=csg_get_property("cg.non-bonded.min");
my $max=csg_get_property("cg.non-bonded.max");
my $step=csg_get_property("cg.non-bonded.step");

my $ndim=$param_N+1;

# Create table with two columns: @r (from grid) and @dummy (0)
my @r;
my @dummy;
my @flag;

open(TMP, ">$tmp");
print TMP "$min 0\n$max 0";
close(TMP);

my @args=("bash","-c","csg_resample --in $tmp --out grid --grid $min:$step:$max");
system(@args);
(readin_table("grid",@r,@dummy,@flag)) || die "$progname: error at readin_table\n";
my @args2=("bash","-c","rm $tmp grid");
system(@args2);

# Read in current simplex table
my (%hash)=readin_simplex_table($infile,$ndim) or die "$progname: error at readin_simplex_table\n";

# Get current parameters
my $sig=${$hash{p_1}}[$p_line_nr];
my $eps=${$hash{p_2}}[$p_line_nr];

# Calculate potential
my @pot;
for (my $i=0;$i<=$#r;$i++){
    # Avoid undefined potential at r=0
    if ($r[$i]>1e-10) {
        $pot[$i]=4*$eps*(($sig/$r[$i])**12-($sig/$r[$i])**6);
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

# Shift potential to zero at the cutoff
for (my $i=0;$i<=$i_cut;$i++){
   $pot[$i]-=$pot[$i_cut];
}

# Save to potential table
saveto_table($outfile,@r,@pot,@flag) or die "$progname: error at saveto_table\n";