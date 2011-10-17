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
#

use strict;
( my $progname = $0 ) =~ s#^.*/##;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
$progname, version %version%
This script:
- calculates the target function (ftar) by comparing the current and target rdf
- sorts ftars according to increasing order of magnitude in simplex table
- does not update if one of the both rdf are undefined

Usage: $progname infile param_N a_line_nr

NEEDS: name kBT

USES: readin_table csg_get_property saveto_table csg_get_property
EOF
  exit 0;
}

die "4 parameters are necessary\n" if ($#ARGV<3);

use CsgFunctions;
use SimplexFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
my $param_N="$ARGV[2]";
my $a_line_nr="$ARGV[3]";

my $name=csg_get_property("cg.non-bonded.name");
my $property=csg_get_property("cg.inverse.simplex.property");
my $weights=csg_get_property("cg.inverse.simplex.weights");

my @property=split(" ", $property);
my @weights=split(" ", $weights);

my $ndim=$param_N+1;

my %hash;
my $ftar_cur;
my $ftar_new=0;
my $p;
my $sum;
# Read in temporary simplex table
foreach $p (0 .. $#property) {
(%hash)=readin_simplex_table("simplex_$name\_$property[$p].tmp",$ndim) or die "$progname: error at readin_simplex_table\n";
$ftar_cur=${$hash{p_0}}[0];
$sum+=$weights[$p]*$ftar_cur;
}
$ftar_new=$sum;

my @ftar;
my @flag;
# Read in temporary simplex table
(%hash)=readin_simplex_table($infile,$ndim) or die "$progname: error at readin_simplex_table\n";

# Define table columns
@ftar=@{$hash{p_0}};
@flag=@{$hash{"p_$ndim"}};

$ftar[$a_line_nr]=$ftar_new;
$flag[$a_line_nr]="complete";

my $mdim=$#ftar+1;

# Save to new simplex table
saveto_simplex_table($outfile,$mdim,$param_N,@ftar,%hash,@flag) or die "$progname: error at saveto_simplex_table\n";
