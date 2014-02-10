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
my $usage="Usage: $progname [OPTIONS] <in> <derivatives_in> <out>";
my $type="non-bonded";
my $sim_prog="none";

while ((defined ($ARGV[0])) and ($ARGV[0] =~ /^-./))
{
  if (($ARGV[0] !~ /^--/) and (length($ARGV[0])>2)){
    $_=shift(@ARGV);
    #short opt having agruments examples fo
    if ( $_ =~ /^-[fo]/ ) {
      unshift(@ARGV,substr($_,0,2),substr($_,2));
    }
    else{
      unshift(@ARGV,substr($_,0,2),"-".substr($_,2));
    }
  }
  if (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help")){
    print <<EOF;
$progname, version %version%
This script converts csg potential files to the tab format
(as read by espresso or lammps or dlpoly).

In addition, it does some magic tricks:
- shift the potential, so that it is zero at the cutoff

$usage

Allowed options:
-h, --help            show this help message
--type XXX            change the type of xvg table
                      Default: $type
--header XXX          Write a special simulation programm header

Examples:
* $progname --type non-bonded table.in table_b0.xvg
EOF
    exit 0;
  }
  elsif ($ARGV[0] eq "--type"){
      shift(@ARGV);
      $type = shift(@ARGV);
  }
  elsif ($ARGV[0] eq "--header"){
      shift(@ARGV);
      $sim_prog = shift(@ARGV);
  }
  else{
    die "Unknown option '".$ARGV[0]."' !\n";
  }
}

#die "$progname: conversion of bonded interaction to generic tables is only implemented for dlpoly!\n" unless (($type eq "non-bonded")||($sim_prog eq "dlpoly"));
die "$progname: conversion of ${type} interaction to generic tables is only implemented for dlpoly!\n" unless (($sim_prog eq "dlpoly"));

die "3 parameters are necessary\n" if ($#ARGV<2);

use CsgFunctions;

my $in_pot="$ARGV[0]";
my $in_deriv_pot="$ARGV[1]";
my $outfile="$ARGV[2]";

my @r;
my @r_repeat;
my @pot;
my @pot_deriv;
my @flag;
my @flag_repeat;
#cutoff is last point
(readin_table($in_pot,@r,@pot,@flag)) || die "$progname: error at readin_table\n";
(readin_table($in_deriv_pot,@r_repeat,@pot_deriv,@flag_repeat)) || die "$progname: error at readin_table\n";

#shift potential so that it is zero at cutoff
for (my $i=0;$i<=$#r;$i++){
   $pot[$i]-=$pot[$#r];
}

open(OUTFILE,"> $outfile") or die "saveto_table: could not open $outfile\n";
# espresso specific header - no other starting comments
if ($sim_prog eq "espresso") {
  printf(OUTFILE "#%d %f %f\n", $#r+1, $r[0],$r[$#r]);
  for(my $i=0;$i<=$#r;$i++){
    printf(OUTFILE "%15.10e %15.10e %15.10e\n",$r[$i],($r[$i]>0)?-$pot_deriv[$i]/$r[$i]:-$pot_deriv[$i], $pot[$i]);
  }
} elsif ($sim_prog eq "lammps") {
  printf(OUTFILE "VOTCA\n");
  printf(OUTFILE "N %i R %f %f\n\n",$#r+1,$r[0],$r[$#r]);
  for(my $i=0;$i<=$#r;$i++){
    printf(OUTFILE "%i %15.10e %15.10e %15.10e\n",$i+1,$r[$i], $pot[$i], -$pot_deriv[$i]);
  }
} elsif ($sim_prog eq "dlpoly") {
  if ($type eq "non-bonded"){
    # see dlpoly manual ngrid = cut/delta + 4  = $#r + 4 as table starts with delta (not 0)
    # number of lines int((ngrid+3)/4)
    for(my $i=0;$i<4*int(($#r+7)/4);$i++){
      printf(OUTFILE "%15.7e",($i>$#r)?0:$pot[$i]);
      printf(OUTFILE "%s",($i%4==3)?"\n":" ");
    }
    for(my $i=0;$i<4*int(($#r+7)/4);$i++){
      # no scaling factor needed 1 kJ/nm *nm = 1*kJ/Angs*Angs
      printf(OUTFILE "%15.7e",($i>$#r)?0:-$pot_deriv[$i]*$r[$i]);
      printf(OUTFILE "%s",($i%4==3)?"\n":" ");
    }
#    printf(OUTFILE "\n");
  } elsif ( $type eq "bond" ) {
#    die "$progname: dlpoly ${type} table not implemented\n";
    for(my $i=0;$i<=$#r;$i++){
      printf(OUTFILE "%15.7e %15.7e %15.7e\n",$r[$i]*10.0, $pot[$i], -$pot_deriv[$i]*$r[$i]);
    }
  } elsif ( $type eq "angle" ||  $type eq "dihedral" ) {
#    die "$progname: dlpoly ${type} table not implemented\n";
    my $RadToDegree=180.0/3.14159265359;
    for(my $i=0;$i<=$#r;$i++){
      printf(OUTFILE "%15.7e %15.7e %15.7e\n",$r[$i]*$RadToDegree, $pot[$i], -$pot_deriv[$i]/$RadToDegree);
    }
  } else {
    #should never happen
    die "$progname: tabulated potentials/forces for dlpoly $type not implemented\n";
  }
} elsif ($sim_prog eq "gromacs") {
  printf(OUTFILE "#This is just a failback, for using different columns use table_to_xvg.pl instead!\n");
  for(my $i=0;$i<=$#r;$i++){
    printf(OUTFILE "%15.10e   %15.10e %15.10e   %15.10e %15.10e   %15.10e %15.10e\n",$r[$i], ,0,0,0,0,$pot[$i], -$pot_deriv[$i]);
  }
} else {
  for(my $i=0;$i<=$#r;$i++){
    printf(OUTFILE "%15.10e %15.10e %15.10e\n",$r[$i], $pot[$i], -$pot_deriv[$i]);
  }
}
close(OUTFILE) or die "Error at closing $outfile\n";

