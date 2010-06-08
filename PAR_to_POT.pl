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
This script flags the current parameter set as 'active'.

Usage: $progname infile outfile param_N p_line_nr

USES: readin_simplex_table saveto_simplex_table csg_get_property csg_resample

NEEDS:

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
my (%hash)=readin_simplex_table($infile,$ndim) or die "$progname: error at readin_simplex_table\n";

# Define columns for ftar and flag
my @ftar=@{$hash{p_0}};
my @flag=@{$hash{"p_$ndim"}};

# Flag current parameter set as 'active'
$flag[$p_line_nr]="active";

my $mdim=$#ftar+1;

# Save to new simplex table
saveto_simplex_table($outfile,$mdim,$param_N,@ftar,%hash,@flag) or die "$progname: error at saveto_simplex_table\n";

my $name=csg_get_property("cg.non-bonded.name");
my $min=csg_get_property("cg.non-bonded.min");
my $max=csg_get_property("cg.non-bonded.max");
my $step=csg_get_property("cg.non-bonded.step");
my $function=csg_get_property("cg.non-bonded.inverse.simplex.function");

# Create table with two columns: @r (from settings) and @dummy (0)
my @r;
my @dummy;

my $tmp=`mktemp tmp_XXX`;
my $grid=`mktemp grid_XXX`;
chop($tmp);
chop($grid);

open(TMP, "> $tmp");
print TMP "$min 0\n$max 0";
close(TMP);

my @args = ("bash", "-c", "csg_resample --in $tmp --out $grid --grid $min:$step:$max");
system(@args);

my $param_string="";
foreach (1..$param_N) {
  $param_string="$param_string ${$hash{\"p_$_\"}}[$p_line_nr]";
}

# Calculate potential
@args = ("bash", "-c", "do_external pot $function $grid $name.pot.new $param_string");
system(@args);