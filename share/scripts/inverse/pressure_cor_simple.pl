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
This script calls the pressure corrections ''\$dU=A*(1-r/r_c)\$'', where
''\$A=-0.1k_B T * \\\\max(1,|p_cur-p_target|*scale) * \\\\sign(p_cur-p_target)\$''

Usage: $progname p_cur outfile kBT min:step:max scale p_target
EOF
  exit 0;
}

die "6 parameters are necessary\n" if ($#ARGV<5);

use CsgFunctions;

my $kBT=$ARGV[2];;
my @range=split(/:/,$ARGV[3]);
defined($range[2]) || die "Not enough number in range $ARGV[3], got ".($#range+1)." need 3\n";
my $max=$range[2];
my $min=$range[0];
my $delta_r=$range[1];
my $scale_factor=$ARGV[4];
my $p_target=$ARGV[5];
my $p_now=$ARGV[0];

#Determine the sign
my $pref;
if ($p_now>$p_target){
   $pref=-0.1*$kBT;
} else {
   $pref=0.1*$kBT;
}

#Determine pressure factor
my $p_factor=($p_now-$p_target)*$scale_factor;
$p_factor=-$p_factor if $p_factor<0;

#Only use pressure factor if not too big
#max is 0.1kbT
$pref*=$p_factor if $p_factor<1;

my @r;
my @pot;
my @flag;
my $outfile="$ARGV[1]";
my $comment="#$progname: p_now=$p_now, p_target=$p_target, prefactor=$pref\n";
for(my $i=$min/$delta_r;$i<=$max/$delta_r;$i++){
  $r[$i]=$i*$delta_r;
  $pot[$i]=$pref*(1-$r[$i]/$max);
  $flag[$i]="i";
}
saveto_table($outfile,@r,@pot,@flag,$comment) || die "$progname: error at save table\n";

