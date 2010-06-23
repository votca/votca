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

die "2 parameters are necessary\n" if ($#ARGV<1);

use CsgFunctions;
use SimplexFunctions;

my $outfile="$ARGV[0]";
my $param_N="$ARGV[1]";

my $name=csg_get_property("cg.non-bonded.name");
my $property=csg_get_property("cg.inverse.simplex.property");
my $weights=csg_get_property("cg.inverse.simplex.weights");

my @property=split(" ", $property);
my @weights=split(" ", $weights);

my $ndim=$param_N+1;

my $count;
my %hash;
my @ftar_cur;
my @ftar_new=();
my @flag_cur;
# Read in temporary simplex table
foreach $count (0 .. $#property) {
(%hash)=readin_simplex_table("simplex_$name\_$property[$count].tmp",$ndim) or die "$progname: error at readin_simplex_table\n";
@ftar_cur=@{$hash{p_0}};
@flag_cur=@{$hash{"p_$ndim"}};
  for (my $i=0;$i<$ndim;$i++) {
    $ftar_new[$i][$count]=$weights[$count]*$ftar_cur[$i];
  }
}

# ------------------- DEFINE TARGET FUNCTION HERE ------------------

my $mdim=$#ftar_new+1;

# Save to new simplex table
saveto_simplex_table($outfile,$mdim,$param_N,@ftar_new,%hash,@flag_cur) or die "$progname: error at saveto_simplex_table\n";
