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
This script marks the current parameter set as active.

Usage: $progname infile outfile param_N p_line_nr

USES: readin_simplex_table saveto_simplex_table

NEEDS: -

EOF
  exit 0;
}

die "3 parameters are nessary\n" if ($#ARGV<2);

use CsgFunctions;
use SimplexFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
my $param_N="$ARGV[2]";
my $p_line_nr="$ARGV[3]";

my $ndim=$param_N+1;

# Read in current simplex table
my (%hash)=readin_simplex_table("$infile","$ndim") or die "$progname: error at readin_simplex_table\n";

# Define first and last column
my @ftar=@{$hash{p_0}};
my @flag=@{$hash{"p_$ndim"}};

# Mark current parameter set as active
$flag[$p_line_nr]="active";

# Save to new simplex table
saveto_simplex_table("$outfile",$param_N,@ftar,%hash,@flag) or die "$progname: error at saveto_simplex_table\n";