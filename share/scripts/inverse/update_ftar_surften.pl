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
This script calculates ftar out of the two rdfs for the 
Simplex Method and arranges the output table in order 
of increasing magnitude of ftar. 
In addtion it does some magic tricks:
- do not update if one of the both rdf are undefined

Usage: $progname target_rdf cur_rdf cur_simplex outfile

NEEDS: cg.inverse.kBT

USES: readin_table saveto_table csg_get_property
EOF
  exit 0;
}

die "5 parameters are necessary\n" if ($#ARGV<4);

use CsgFunctions;
use SimplexFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
my $param_N="$ARGV[2]";
my $a_line_nr="$ARGV[3]";
my $prop_N="$ARGV[4]";

my $property="surften";
my $name=csg_get_property("cg.non-bonded.name");
my $sim_prog=csg_get_property("cg.inverse.program");

# Get tgt surften
my $surften_tgt=csg_get_interaction_property("inverse.simplex.surften.target");

# Calculate new surften
my @args = ("bash", "-c", "for_all non-bonded do_external $property $sim_prog");
system(@args);
undef(@args);

my $surften_cur;
open(SURFTEN_CUR, "<$name.surften.cur");
while (<SURFTEN_CUR>) {
  $surften_cur=$_;
}
close(SURFTEN_CUR);

my @ftar_cur;
my @flag_cur;

my $ndim=$param_N+1;

# Read in temporary simplex table
my (%hash)=readin_simplex_table($infile,$ndim) or die "$progname: error at readin_simplex_table\n";

# Define table columns
@ftar_cur=@{$hash{p_0}};
@flag_cur=@{$hash{"p_$ndim"}};

# ------------------- DEFINE TARGET FUNCTION HERE ------------------
# Calculate ftar
my $ftar=abs(($surften_cur-$surften_tgt)/$surften_tgt)*100;

my $mdim;
if ($prop_N == 1) {
  $ftar_cur[$a_line_nr]=$ftar;
  $flag_cur[$a_line_nr]="complete";
  $mdim=$#ftar_cur+1;
}
else {
  $ftar_cur[0]=$ftar;
  for (my $j=1;$j<=$param_N;$j++){
    ${$hash{"p_$j"}}[0]=${$hash{"p_$j"}}[$a_line_nr];
  }
  $flag_cur[0]="complete";
  $mdim=1;
}

# Save to new simplex table
saveto_simplex_table($outfile,$mdim,$param_N,@ftar_cur,%hash,@flag_cur) or die "$progname: error at saveto_simplex_table\n";
