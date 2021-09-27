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
my $usage="Usage: $progname [OPTIONS] kbint target_kbint outfile kBT min:step:max int_start:int_end ramp_factor";

while ((defined ($ARGV[0])) and ($ARGV[0] =~ /^-./))
{
  if (($ARGV[0] !~ /^--/) and (length($ARGV[0])>2)){
    $_=shift(@ARGV);
    #short opt having agruments examples fo
    if ( $_ =~ /^-[fo]/ ) {
      unshift(@ARGV,substr($_,0,2),substr($_,2));
    } else{
      unshift(@ARGV,substr($_,0,2),"-".substr($_,2));
    }
  }
  if (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help")){
     print <<END;
$progname, version %version%

This script calculates Kirkwood-Buff correction as described in:
P. Ganguly, D. Mukherji, C. Junghans, N. F. A. van der Vegt,
Kirkwood-Buff coarse-grained force fields for aqueous solutions,
J. Chem. Theo. Comp., 8, 1802 (2012), doi:10.1021/ct3000958

$usage

Allowed options:
-h, --help            Show this help message
END
    exit 0;
  }else{
    die "Unknown option '".$ARGV[0]."' !\n";
  }
}

die "7 parameters are necessary\n" if ($#ARGV<6);

use CsgFunctions;

my $kbt=$ARGV[3];
my @irange=split(/:/,$ARGV[5]);
defined($irange[1]) || die "Not enough number in irange $ARGV[5], got ".($#irange+1)." need 2\n";
my $int_start=$irange[0];
my $int_stop=$irange[1];
my $ramp_factor=$ARGV[6];

my @range=split(/:/,$ARGV[4]);
defined($range[2]) || die "Not enough number in range $ARGV[4], got ".($#range+1)." need 3\n";
my $r_min=$range[0];
my $r_ramp=$range[2];
my $delta_r=$range[1];

my $aim_kbint_file="$ARGV[0]";
my @r_aim;
my @kbint_aim;
my @flags_aim;
(readin_table($aim_kbint_file,@r_aim,@kbint_aim,@flags_aim)) || die "$progname: error at readin_table\n";

my $cur_kbint_file="$ARGV[1]";
my @r_cur;
my @kbint_cur;
my @flags_cur;
(readin_table($cur_kbint_file,@r_cur,@kbint_cur,@flags_cur)) || die "$progname: error at readin_table\n";

#should never happen due to resample, but better check
die "Different grids \n" if (($r_aim[1]-$r_aim[0]-$r_cur[1]+$r_cur[0])>0.0001);
die "Different start potential point \n" if (($r_aim[0]-$r_cur[0]) > 0.0001);
die "Different end potential point \n" if ( $#r_aim != $#r_cur );

my $j=0;
my $avg_int=0;
for (my $i=0;$i<=$#r_aim;$i++){
  if (($r_aim[$i]>=$int_start) && ($r_aim[$i]<=$int_stop)) {
     $avg_int+=$kbint_cur[$i]-$kbint_aim[$i];
     $j++;
  }
}
$avg_int/=$j;

my $comment="#$progname: avg_int($int_start:$int_stop)=$avg_int ramp_factor=$ramp_factor r_ramp=$r_ramp\n";
my @dpot;
my @flag;
for (my $i=0;$i<=$#r_aim;$i++){
  if ($r_aim[$i]> $r_ramp) {
    $dpot[$i]=0; #beyond r_ramp correction is 0
  } else {
    $dpot[$i]=($avg_int*$ramp_factor*(1.0-($r_aim[$i]/$r_ramp)))*$kbt;
  }
  $flag[$i]="i";
}

my $outfile="$ARGV[2]";
saveto_table($outfile,@r_aim,@dpot,@flag,$comment) || die "$progname: error at save table\n";
