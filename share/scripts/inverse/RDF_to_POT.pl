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
$progname, %version%
This script converts rdf to pot of mean force (F(r)=-k_B T*ln(g(r))

In addtion, it does some magic tricks:
- do not crash when calc log(0)
- extrapolate the beginning of pot
- the maximum to interpolate is pot_max (see xml)
- bigger value will be set to that max
- shift the potential, so that it is zero at the cutoff
- set all values to zero after the cutoff

Usage: $progname infile outfile

USES: readin_table csg_get_property csg_get_property csg_get_interaction_property saveto_table

NEEDS: cg.inverse.kBT max

EOF
  exit 0;
}

die "2 parameters are nessary\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my $pref=csg_get_property("cg.inverse.kBT");
my $r_cut=csg_get_interaction_property("max");

my @r;
my @rdf;
my @flag;
(readin_table($infile,@r,@rdf,@flag)) || die "$progname: error at readin_table\n";

my @pot;
for (my $i=0;$i<=$#r;$i++){
#  if ($flag[$i] eq "i"){
    #rdf = 0 will give undefined pot
    if ($rdf[$i]>1e-10) {
      $pot[$i]=-$pref*log($rdf[$i]);
    }
    else {
      $pot[$i]="nan";
      $flag[$i]="u";
    }
#  }
}

#find first defined value (begining for r=0)
#but it is more stable to search first undefined value begin
#beginning form large r
my $first_undef_bin=-1;
for (my $i=$#pot;$i>=0;$i--){
   if ($flag[$i] eq "u") {
     $first_undef_bin=$i;
     last;
   }
}

#find i which is the cutoff
my $i_cut=$#r;
for (my $nr=0;$nr<=$#r;$nr++){
   if ($r[$nr]>=$r_cut) {
     $i_cut=$nr;
     last;
   }
}

#shift potential so that it is zero at cutoff
#first do the shift, then extrapolation
for (my $i=0;$i<=$i_cut;$i++){
   $pot[$i]-=$pot[$i_cut] unless  ($flag[$i] =~ /[u]/);
}

#quadratic extrapolation at the begining
my $slope=$pot[$first_undef_bin+1]-$pot[$first_undef_bin+2];
for (my $i=$first_undef_bin;$i>=0;$i--){
   $slope+=$slope;
   $pot[$i]=$pot[$i+1]+$slope;
   $flag[$i]="o";
}

# set end of the potential to zero
for (my $i=$i_cut;$i<$#flag;$i++) {
  $pot[$i]=0;
  $flag[$i]="o";
}

saveto_table($outfile,@r,@pot,@flag) || die "$progname: error at save table\n";

