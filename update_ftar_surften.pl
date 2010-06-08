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

die "3 parameters are necessary\n" if ($#ARGV<2);

use CsgFunctions;
use SimplexFunctions;

my $name=csg_get_property("cg.non-bonded.name");
my $surften_tgt=csg_get_interaction_property("inverse.target");

my $surften_cur;
open(SURFTEN_CUR, "<surften.cur");
while (<SURFTEN_CUR>) {
  $surften_cur=$_;
}
close(SURFTEN_CUR);

# Create an empty rdf file
open (RDF, "> $name.dist.new") || die "Could not open file $_[0]\n";
close(RDF);

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
my $param_N="$ARGV[2]";
my $a_line_nr="$ARGV[3]";

my @ftar_cur;
my @sig_cur;
my @eps_cur;
my @flag_cur;

my $ndim=$param_N+1;

# Read in temporary simplex table
my (%hash)=readin_simplex_table($infile,$ndim) or die "$progname: error at readin_simplex_table\n";

# Define table columns
@ftar_cur=@{$hash{p_0}};
@flag_cur=@{$hash{"p_$ndim"}};

# --------------------- DEFINE PARAMETERS HERE ---------------------
@sig_cur=@{$hash{p_1}};
@eps_cur=@{$hash{p_2}};

# ------------------- DEFINE TARGET FUNCTION HERE ------------------
# Calculate ftar

$ftar_cur[$a_line_nr]=abs(($surften_cur-$surften_tgt)/$surften_tgt);

my @args=("bash","-c","echo $ftar_cur[$a_line_nr]");
system(@args);

my @ftar_new;
@ftar_new=@ftar_cur;

# Flag current parameter set as 'complete'
$flag_cur[$a_line_nr]="complete";

my $mdim=$#ftar_cur+1;

# Save to new simplex table
saveto_simplex_table($outfile,$mdim,$param_N,@ftar_cur,%hash,@flag_cur) or die "$progname: error at saveto_simplex_table\n";
