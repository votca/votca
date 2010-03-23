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

my $aim_rdf_file="$ARGV[0]";
my @r_aim;
my @rdf_aim;
my @flags_aim;
(readin_table($aim_rdf_file,@r_aim,@rdf_aim,@flags_aim)) || die "$progname: error at readin_table\n";

my $cur_rdf_file="$ARGV[1]";
my @r_cur;
my @rdf_cur;
my @flags_cur;
(readin_table($cur_rdf_file,@r_cur,@rdf_cur,@flags_cur)) || die "$progname: error at readin_table\n";

my $cur_ftar_file="$ARGV[2]";
my @ftar_cur;
my @sig_cur;
my @eps_cur;
my @flag_simplex;
(readin_simplex_table($cur_ftar_file,@ftar_cur,@sig_cur,@eps_cur,@flag_simplex)) || die "$progname: error at readin_simplex_table\n";

my $new_ftar_file="$ARGV[3]";
my $a_line_nr="$ARGV[4]";

# Should never happen due to resample, but better check
die "Different grids \n" if (($r_aim[1]-$r_aim[0])!=($r_cur[1]-$r_cur[0]));
die "Different start point \n" if (($r_aim[0]-$r_cur[0]) > 0.0);

$ftar_cur[$a_line_nr]=calc_ftar($a_line_nr,@r_cur,@rdf_aim,@rdf_cur,@sig_cur,@eps_cur);

my @ftar_new;
@ftar_new=@ftar_cur;

$ftar_new[$a_line_nr]="complete";

saveto_table($new_ftar_file,@ftar_new,@sig_cur,@eps_cur) || die "$progname: error at save table\n";