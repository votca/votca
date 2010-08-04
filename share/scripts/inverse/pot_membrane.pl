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

Usage: $progname infile outfile param_cur

USES: readin_table saveto_table

NEEDS:
EOF
  exit 0;
}

# Give number of parameters
if (defined($ARGV[0])&&("$ARGV[0]" eq "--nparams")){
  print "2\n";
  exit 0;
}

die "3 parameters are nessary\n" if ($#ARGV<2);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
# --------------------- DEFINE PARAMETERS HERE ---------------------
my $sig="$ARGV[2]";
my $eps="$ARGV[3]";

# Read in empty potential table
my @r;
my $r_cut_rep=(2**(1/6))*$sig;
my $r_cut_tot=csg_get_interaction_property("max");
my $w_cut=$r_cut_tot-$r_cut_rep;
my $pi=3.14159265;
my @pot;
my @pot_rep;
my @pot_att;
my @flag;
(readin_table($infile,@r,@pot,@flag)) || die "$progname: error at readin_table\n";

# -------------------- DEFINE POTENTIAL HERE -----------------------
# Calculate potential
for (my $i=0;$i<=$#r;$i++){
    # Avoid undefined potential at r=0
    if ($r[$i]>1e-10) {

      # Attractive part
      if ($r[$i]<$r_cut_rep) {
        $pot_att[$i]=-$eps;
      }
      else {
        $pot_att[$i]=-$eps*(cos($pi*($r[$i]-$r_cut_rep)/(2*$w_cut)))**2;        
      }
 
      # Repulsive part
      if ($r[$i]<=$r_cut_rep) {
        $pot_rep[$i]=(4*$eps*((($sig/$r[$i])**12)-(($sig/$r[$i])**6)+(1/4)));
      }
      else {
        $pot_rep[$i]=0;
      }

    # Total potential
    $pot[$i]=$pot_rep[$i]+$pot_att[$i];
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
