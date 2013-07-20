#! /usr/bin/perl -w
#
# Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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
This script converts csg potential files to the tab format
(as read by espresso and lammps).

In addition, it does some magic tricks:
- shift the potential, so that it is zero at the cutoff

Usage: $progname in_pot in_deriv_pot outfile
EOF
  exit 0;
}

die "3 parameters are necessary\n" if ($#ARGV<2);

use CsgFunctions;

my $in_pot="$ARGV[0]";
my $in_deriv_pot="$ARGV[1]";
my $outfile="$ARGV[2]";

my $sim_prog=csg_get_property("cg.inverse.program");

my @r;
my @r_repeat;
my @pot;
my @d_pot;
my @flag;
my @flag_repeat;
#cutoff is last point
(readin_table($in_pot,@r,@pot,@flag)) || die "$progname: error at readin_table\n";
(readin_table($in_deriv_pot,@r_repeat,@d_pot,@flag_repeat)) || die "$progname: error at readin_table\n";

#shift potential so that it is zero at cutoff
for (my $i=0;$i<=$#r;$i++){
   $pot[$i]-=$pot[$#r];
}

my @minus_force=@d_pot;

# Smooth out force (9-point avg) 
for (my $i=4;$i<$#r_repeat-3;$i++){
		$minus_force[$i]=($d_pot[$i-4]+$d_pot[$i-3]+$d_pot[$i-2]
								+$d_pot[$i-1]+$d_pot[$i]+$d_pot[$i+1]+$d_pot[$i+2]
								+$d_pot[$i+3]+$d_pot[$i+4])/(9.);
}

if ($sim_prog eq "espresso") {
  # add extra 1/r factor for ESPResSo
  for (my $i=0;$i<=$#r_repeat;$i++){
		$minus_force[$i]*=1.0/$r_repeat[$i] if ($r_repeat[$i] > 0);
  } 
}

open(OUTFILE,"> $outfile") or die "saveto_table: could not open $outfile\n";
# espresso specific header - no other starting comments
if ($sim_prog eq "espresso") {
  printf(OUTFILE "#%d %f %f\n", $#r+1, $r[0],$r[$#r]);
  for(my $i=0;$i<=$#r;$i++){
    printf(OUTFILE "%15.10e %15.10e %15.10e\n",$r[$i], -$minus_force[$i], $pot[$i]);
  }
} elsif ($sim_prog eq "lammps") {
  printf(OUTFILE "VOTCA\n");
  printf(OUTFILE "N %i R %f %f\n\n",$#r+1,$r[0],$r[$#r]);
  for(my $i=0;$i<=$#r;$i++){
    printf(OUTFILE "%i %15.10e %15.10e %15.10e\n",$i+1,$r[$i], $pot[$i], -$minus_force[$i]);
  }
} else {
  for(my $i=0;$i<=$#r;$i++){
    printf(OUTFILE "%i %15.10e %15.10e %15.10e\n",$r[$i], $pot[$i], -$minus_force[$i]);
  }
}
close(OUTFILE) or die "Error at closing $outfile\n";

