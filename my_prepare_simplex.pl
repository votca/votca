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
This script takes the input parameters matrix and adds
columns for ftar and flag, creating the standard simplex
table file.

Usage: $progname infile outfile statefile param_N

USES: readin_init_simplex_table saveto_simplex_table

NEEDS: -

EOF
  exit 0;
}

die "4 parameters are nessary\n" if ($#ARGV<3);

use SimplexFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
my $statefile="$ARGV[2]";
my $param_N="$ARGV[3]";

# Create a state file
open (STATE, "> $statefile") || die "Could not open file $_[0]\n";
print STATE "Transformation=None\n";
close STATE;

my $ndim=$param_N+1;

# Read in current simplex table
my (%hash)=readin_init_simplex_table($infile,$param_N) or die "$progname: error at readin_simplex_table\n";

my @ftar;
my @flag;

for my $i (0 .. $param_N) {
   $ftar[$i]="0.0";
   $flag[$i]="pending";
}

# Save to new simplex table
saveto_simplex_table($outfile,$param_N,@ftar,%hash,@flag) or die "$progname: error at saveto_simplex_table\n";