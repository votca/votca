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

die "2 parameters are nessary\n" if ($#ARGV<1);

use SimplexFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my @sig;
my @eps;

readin_init_simplex_table($infile,@sig,@eps) || die "$progname: error at readin_init_simplex_table\n";

my @f_target;
my @flag_simplex;

for my $i (0 .. $#sig) {
   $f_target[$i]="0";
   $flag_simplex[$i]="pending";
}

my $p_nr=$#sig;

# Create a state file
open (STATE, ">state.new") || die "Could not open file $_[0]\n";
print STATE "Transformation=None\n";
close STATE;

saveto_simplex_table($outfile,@f_target,@sig,@eps,@flag_simplex) || die "$progname: error at saveto_simplex_table\n";