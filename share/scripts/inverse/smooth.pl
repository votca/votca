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
This script smoothes a table

Usage: $progname infile outfile
EOF
  exit 0;
}

die "2 parameters are nessary\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my @r_cur;
my @pot_cur;
my @flag_cur;
(readin_table($infile,@r_cur,@pot_cur,@flag_cur)) || die "$progname: error at readin_table\n";

my $outfile="$ARGV[1]";
my @pot;

# TODO: think about addition rules
# now I did it like that to always maintain interval of interest in all potentials
for (my $i=1;$i<$#r_cur;$i++){
  $pot[$i]=$pot_cur[$i];
  if($flag_cur[$i] eq "i") {
    $pot[$i] = 0.25*$pot_cur[$i-1] + 0.5*$pot_cur[$i] + 0.25*$pot_cur[$i+1];
  }
}

$pot[0]=$pot_cur[0];
$pot[$#pot_cur]=$pot_cur[$#pot_cur];

if($flag_cur[0] eq "i") {
  $pot[0] = (2.*$pot_cur[0] + $pot_cur[1])/3;
}
if($flag_cur[$#pot_cur] eq "i") {
  $pot[$#pot_cur] = (2.*$pot_cur[$#pot_cur] + $pot_cur[$#pot_cur-1])/3;
}

saveto_table($outfile,@r_cur,@pot,@flag_cur) || die "$progname: error at save table\n";

