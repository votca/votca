#! /usr/bin/perl -w
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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
This script converts rdf to pot of mean force (''\$F(r)=-k_B T\\\\ln g(r)\$'')

In addtion, it does some magic tricks:
- do not crash when calc log(0)

Usage: $progname infile outfile
EOF
  exit 0;
}

die "2 parameters are necessary\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my $pref=csg_get_property("cg.inverse.kBT");
my $rdf_min=csg_get_property("cg.inverse.rdf_min");

my @r;
my @rdf;
my @flag;
(readin_table($infile,@r,@rdf,@flag)) || die "$progname: error at readin_table\n";

my @pot;
for (my $i=0;$i<=$#r;$i++){
    if ($rdf[$i]>$rdf_min) {
      $pot[$i]=-$pref*log($rdf[$i]);
    }
    else {
      $pot[$i]="nan";
      $flag[$i]="u";
    }
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
die "All data points from file '$infile' are invalid after Boltzmann inversion, please check if your distribution is a valid rdf.\n" if ($first_undef_bin==$#pot);

#set point at beginning to invalid
for (my $i=$first_undef_bin;$i>=0;$i--){
   $pot[$i]=$pot[$first_undef_bin+1];
   $flag[$i]="o";
}

saveto_table($outfile,@r,@pot,@flag) || die "$progname: error at save table\n";
