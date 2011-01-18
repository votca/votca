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
This script adds up two potentials
In addition, it does some magic tricks:
- order of infiles MATTER !!!!
- if infile2 contains an undefined value, it uses the value from infile1
- if value for infile1 and infile2 are both invalid, the result is also invalid

Usage: $progname infile1 infile2 outfile

NEEDS:

USES: readin_table saveto_table
EOF
  exit 0;
}

die "3 parameters are nessary\n" if ($#ARGV<2);

use CsgFunctions;

my $infile="$ARGV[0]";
my @r_cur;
my @pot_cur;
my @flag_cur;
(readin_table($infile,@r_cur,@pot_cur,@flag_cur)) || die "$progname: error at readin_table\n";

my $infile2="$ARGV[1]";
#delta is just a name
my @r_delta;
my @pot_delta;
my @flag_delta;
(readin_table($infile2,@r_delta,@pot_delta,@flag_delta)) || die "$progname: error at readin_table\n";

#should never happen, but ....
die "Different grids\n" if (($r_delta[1]-$r_delta[0]-$r_cur[1]+$r_cur[0])>0.0001);
die "Different start point \n" if (($r_delta[0]-$r_cur[0]) > 0.0);

my $outfile="$ARGV[2]";
my @pot;
my @flag;

# TODO: think about addition rules
# now I did it like that to always maintain interval of interest in all potentials
# shount that just be a < instead of <= ??
for (my $i=0;$i<=$#r_cur;$i++){
  if($flag_cur[$i] eq "u" || $flag_delta[$i] eq "u") {
    $pot[$i] = $pot_cur[$i];  # is already nan or we don't change
    $flag[$i] = "u";
  }
  else {
    $pot[$i]=$pot_cur[$i]+$pot_delta[$i];
    $flag[$i] = $flag_cur[$i];
  }
  #if ($flag_cur[$i] eq "i"){
  #  if ($flag_delta[$i] eq "i"){
  #    $pot[$i]=$pot_cur[$i]+$pot_delta[$i];
  #  } else {
  #    $pot[$i]=$pot_cur[$i];
  #  }
  #  $flag[$i]="i";
  #} else {
  #  $pot[$i]="nan";
  #  $flag[$i]="u";
  #}
}
saveto_table($outfile,@r_cur,@pot,@flag) || die "$progname: error at save table\n";

