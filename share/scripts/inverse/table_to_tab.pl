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
This script convert csg potential files to tab format (as read by espresso)
Potential are copy in the C12 column
In addtion it does some magic tricks:
- bigger value will be set to pot_max (see xml)
- shift the potential, so that it is zero at the cutoff
- set all values to zero after the cutoff

Usage: $progname infile outfile

NEEDS: cg.inverse.espresso.pot_max cg.inverse.espresso.table_end cg.inverse.espresso.table_bins

USES: csg_get_property saveto_table readin_table
EOF
  exit 0;
}

die "2 parameters are nessary\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my $espresso_max=csg_get_property("cg.inverse.espresso.pot_max");
my $table_end=csg_get_property("cg.inverse.espresso.table_end");
my $table_bins=csg_get_property("cg.inverse.espresso.table_bins");

my @r;
my @pot;
my @flag;
(readin_table($infile,@r,@pot,@flag)) || die "$progname: error at readin_table\n";

#Avoid very large numbers
for (my $i=0;$i<=$#r;$i++) {
  $pot[$i]=$espresso_max if $pot[$i]>$espresso_max;
  $pot[$i]=-$espresso_max if $pot[$i]<-$espresso_max;
}

#cutoff is last point
my $i_cut=$#r;

#shift potential so that it is zero at cutoff
for (my $i=0;$i<=$i_cut;$i++){
   $pot[$i]-=$pot[$i_cut];
}

# set end of the potential to zero
for (my $i=$i_cut;$i<=$table_end/$table_bins;$i++) {
  $pot[$i]=0;
  $r[$i]=$r[$i-1]+$table_bins;
}

my @force;

#calc force
$force[0]=0;
for (my $i=1;$i<$#r;$i++){
   $force[$i]=-($pot[$i+1]-$pot[$i-1])/($r[$i+1]-$r[$i-1]);
}
$force[$#r]=0.0;

open(OUTFILE,"> $outfile") or die "saveto_table: could not open $outfile\n";
# espresso specific header - no other starting comments
my $num_bins = $table_end/$table_bins;
printf(OUTFILE "#%d 0 %f\n", $num_bins, $table_end);
for(my $i=0;$i<=$#r;$i++){
  printf(OUTFILE "%15.10e %15.10e %15.10e\n",
    $r[$i], $pot[$i], $force[$i]);
}
close(OUTFILE) or die "Error at closing $outfile\n";

