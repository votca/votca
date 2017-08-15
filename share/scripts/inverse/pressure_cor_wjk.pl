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
This script calls the pressure corrections  like in
Wan, Junghans & Kremer, Euro. Phys. J. E 28, 221 (2009)
Basically dU=A*(1-r/r_c) with A= -max(0.1k_B T, Int ) * sign(p_cur-p_target)
and Int is the integral from Eq. 7 in the paper.

Usage: $progname p_cur outfile kBT min:step:max scale p_target particle_dens rdf_file
EOF
  exit 0;
}

die "8 parameters are necessary\n" if ($#ARGV<7);


use CsgFunctions;

my $kBT=$ARGV[2];;
my @range=split(/:/,$ARGV[3]);
defined($range[2]) || die "Not enough number in range $ARGV[3], got ".($#range+1)." need 3\n";
my $max=$range[2];
my $min=$range[0];
my $delta_r=$range[1];

my $partDens=$ARGV[6];
my $scale_factor=$ARGV[4];

my $pi= 3.14159265;
my $bar_to_SI = 0.06022; # 1bar=0.06022 kJ/(nm mol)

my $p_target=$ARGV[5];
my $p_now=$ARGV[0];

# load current rdf
my $cur_rdf_file="$ARGV[7]";
my @r_cur;
my @rdf_cur;
my @flags_cur;

(readin_table($cur_rdf_file,@r_cur,@rdf_cur,@flags_cur)) || die "$progname: error at readin_table\n";

# calculate prefactor from rdf
my $integral=0.0;
my $x;
for(my $i=1;$i<$max/$delta_r;$i++){
	$x=$i*$delta_r;
	$integral+=$x*$x*$x*$delta_r*$rdf_cur[$i];
}
my $pref;

$integral += ($delta_r/2*$rdf_cur[$max/$delta_r]*$max*$max*$max);
$pref = -3*$max*($p_now-$p_target)*$bar_to_SI;
$pref /= 2*$pi*$partDens*$partDens*$integral;

# use max($pref, +-0.1kt) as prefactor

my $temp;
$temp=$pref;
$temp = -1*$temp if $temp<0;
if ($temp > 0.1*$kBT){
	if ($pref >0){
		$pref=0.1*$kBT;
	}else{
		$pref=-0.1*$kBT;
	}
}

$pref=$pref*$scale_factor;
print "Pressure correction factor: A=$pref\n";

# my $prefile="${name}.pressure.prefactor";
# saveto_table($prefile,$pref) || die "$progname: error at save table\n";

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

