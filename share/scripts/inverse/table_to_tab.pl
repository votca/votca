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

Usage: $progname in_pot in_deriv_pot outfile

NEEDS: cg.inverse.espresso.pot_max cg.inverse.espresso.table_end cg.inverse.espresso.table_bins

USES: csg_get_property saveto_table readin_table
EOF
  exit 0;
}

die "3 parameters are nessary\n" if ($#ARGV<2);

use CsgFunctions;

my $in_pot="$ARGV[0]";
my $in_deriv_pot="$ARGV[1]";
my $outfile="$ARGV[2]";

my $table_end=csg_get_property("cg.inverse.espresso.table_end");
my $table_bins=csg_get_property("cg.inverse.espresso.table_bins");

my @r;
my @r_repeat;
my @pot;
my @d_pot;
my @flag;
my @flag_repeat;
(readin_table($in_pot,@r,@pot,@flag)) || die "$progname: error at readin_table\n";
(readin_table($in_deriv_pot,@r_repeat,@d_pot,@flag_repeat)) || die "$progname: error at readin_table\n";

#cutoff is last point
my $i_cut=$#r;

#shift potential so that it is zero at cutoff
for (my $i=0;$i<=$i_cut;$i++){
   $pot[$i]-=$pot[$i_cut];
}

my @force=@d_pot;

# set end of the potential to zero
for (my $i=$i_cut;$i<=$table_end/$table_bins;$i++) {
  $pot[$i]=0;
	$force[$i]=0;
  $r[$i]=$r[$i-1]+$table_bins;
}


# Smooth out force (9-point avg) 
for (my $i=4;$i<$#r_repeat-3;$i++){
		$force[$i]=($d_pot[$i-4]+$d_pot[$i-3]+$d_pot[$i-2]
								+$d_pot[$i-1]+$d_pot[$i]+$d_pot[$i+1]+$d_pot[$i+2]
								+$d_pot[$i+3]+$d_pot[$i+4])/(9.);
}
# add extra 1/r factor for ESPResSo
for (my $i=1;$i<$#r_repeat;$i++){
		$force[$i]*=-1.0/$r_repeat[$i];
}
$force[0]=$force[1];

open(OUTFILE,"> $outfile") or die "saveto_table: could not open $outfile\n";
# espresso specific header - no other starting comments
my $num_bins = $table_end/$table_bins;
printf(OUTFILE "#%d 0 %f\n", $num_bins, $table_end);
for(my $i=0;$i<=$#r;$i++){
  printf(OUTFILE "%15.10e %15.10e %15.10e\n",
    $r[$i], $force[$i], $pot[$i]);
}
close(OUTFILE) or die "Error at closing $outfile\n";

