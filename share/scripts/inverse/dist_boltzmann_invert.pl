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
my $usage="Usage: $progname [OPTIONS] <in> <out>";

my $kbT=undef;
my $rdf_min=1e-10;
my $type="non-bonded";

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
  if (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help"))
  {
    print <<END;
$progname, version %version%
Boltzmann inverts a distribution (''\$F(r)=-k_B T\\\\ln g(r)\$'')

In addtion, it does some magic tricks:
- do not crash when calc log(0)
- choose the right normalization depending on the type of interaction

$usage

Allowed options:
    --kbT NUMBER      use NUMBER as ''\$k_B*T\$'' for the entropic part
    --type XXX        change the type of interaction
                      Default: $type
    --min XXX         minimum value to consider
                      Default: $rdf_min
-h, --help            Show this help message

Possible types: non-bonded, bond, angle, dihedral, bonded

Examples:
* $progname --kbT 2.49435 --min 0.001 tmp.dist tmp.pot
END
    exit;
  }
  elsif ($ARGV[0] eq "--kbT"){
    shift(@ARGV);
    $kbT=shift(@ARGV);
  }
  elsif ($ARGV[0] eq "--type"){
    shift(@ARGV);
    $type = shift(@ARGV);
  }
  elsif ($ARGV[0] eq "--min"){
    shift(@ARGV);
    $rdf_min = shift(@ARGV);
  }
  else {
    die "Unknown option '".$ARGV[0]."' !\n";
  }
}
die "2 parameters are necessary\n" if ($#ARGV<1);

die "$progname: kbT not defined specify it with --kbT option\n" unless defined($kbT);
die "$progname: conversion of bonded interaction to generic tables is not implemented yet!" unless ($type eq "non-bonded");

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my @r;
my @rdf;
my @flag;
(readin_table($infile,@r,@rdf,@flag)) || die "$progname: error at readin_table\n";

my @pot;
for (my $i=0;$i<=$#r;$i++){
    if ($rdf[$i]>$rdf_min) {
      $pot[$i]=-$kbT*log($rdf[$i]);
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
