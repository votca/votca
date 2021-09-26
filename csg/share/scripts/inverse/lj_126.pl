#! /usr/bin/perl -w
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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
This script calculates the LJ 12-6 potential ''\$U=C12/r^12 - C6/r^6\$''

Usage: $progname outfile min:step:max C6 C12
EOF
  exit 0;
}

die "4 parameter is necessary\n" if ($#ARGV<3);

use CsgFunctions;

my @range=split(/:/,$ARGV[1]);
defined($range[2]) || die "Not enough number in range $ARGV[1], got ".($#range+1)." need 3\n";
my $max=$range[2];
my $min=$range[0];
my $delta_r=$range[1];

my $c12=$ARGV[3];
my $c6=$ARGV[2];

my @r;
my @pot;
my @flag;
my $outfile="$ARGV[0]";
my $comment="#$progname: C12=$c12, C6=$c6\n";
for(my $i=$min/$delta_r;$i<=$max/$delta_r;$i++){
  $r[$i]=$i*$delta_r;
  if($r[$i]>0.0){
      $pot[$i]=$c12/(($r[$i])**12) - $c6/(($r[$i])**6);
  } else {
      $pot[$i]=1.0E20; # very large number
  }
  $flag[$i]="i";
}
saveto_table($outfile,@r,@pot,@flag,$comment) || die "$progname: error at save table\n";
