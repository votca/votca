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
This script calculates the jackknife error from existing tables
* full     = table calculated with full dataset
* blocks   = tables calculated with 1 block missing
* outfile  = file to write results

Usage: $progname out full block1 block2 ...

USES: readin_table saveto_table_err

EOF
  exit 0;
}

die "3 parameters are nessary\n" if ($#ARGV<2);

use CsgFunctions;

my $file_full="$ARGV[1]";
my @r_full;
my @val_full;
my @flag_full;

my $outfile="$ARGV[0]";
my @err;


(readin_table($file_full,@r_full,@val_full,@flag_full)) || die "$progname: error at readin_table\n";

for (my $i=0;$i<=$#r_full;$i++) {
  $err[$i]=0;
}

shift @ARGV;
shift @ARGV;

my $nblocks = 0;
while (@ARGV > 0) {
  my $file_cur="$ARGV[0]";
  my @r_cur;
  my @val_cur;
  my @flag_cur;

  (readin_table($file_cur,@r_cur,@val_cur,@flag_cur)) || die "$progname: error at readin_table\n";
  #should never happen, but ....
  #die "Different grids\n" if (($r_delta[1]-$r_delta[0]-$r_cur[1]+$r_cur[0])>0.0001);
  #die "Different start point \n" if (($r_delta[0]-$r_cur[0]) > 0.0);

  for (my $i=0;$i<=$#r_cur;$i++) {
      $err[$i] += ($val_cur[$i] - $val_full[$i])**2;  # is already nan or we don't change
  }
  shift @ARGV;
  $nblocks = $nblocks + 1;
}

for (my $i=0;$i<=$#r_full;$i++) {
  $err[$i]=sqrt(($nblocks-1)/$nblocks*$err[$i]);
}

saveto_table_err($outfile,@r_full,@val_full,@flag_full,@err) || die "$progname: error at save table\n";

