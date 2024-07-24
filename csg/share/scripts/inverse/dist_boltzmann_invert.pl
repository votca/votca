#! /usr/bin/env -S perl -w
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
my $dist_min=1e-10;
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
Boltzmann inverts a distribution (''\$F(x)=-k_B T\\\\ln g(x)\$'')

In addtion, it does some magic tricks:
- do not crash when calc log(0)
- choose the right normalization depending on the type of interaction
- input dist should be unnormalized (like csg_stat calcs it)

$usage

Allowed options:
    --kbT NUMBER      use NUMBER as ''\$k_B*T\$'' for the entropic part
    --type XXX        change the type of interaction
                      Default: $type
    --min XXX         minimum value to consider
                      Default: $dist_min
-h, --help            Show this help message

Possible types: non-bonded, bond, angle, dihedral

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
    die "$progname: Unsupported type of interatction: $type -> go and implement it\n" unless (( "$type" eq "bond" ) or ("$type" eq "dihedral") or ("$type" eq "angle") or ("$type" eq "non-bonded"));
  }
  elsif ($ARGV[0] eq "--min"){
    shift(@ARGV);
    $dist_min = shift(@ARGV);
  }
  else {
    die "Unknown option '".$ARGV[0]."' !\n";
  }
}
die "2 parameters are necessary\n" if ($#ARGV<1);

die "$progname: kbT not defined specify it with --kbT option\n" unless defined($kbT);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my @x;
my @dist;
my @flag;
(readin_table($infile,@x,@dist,@flag)) || die "$progname: error at readin_table\n";

my @pot;
for (my $i=0;$i<=$#x;$i++){
    if ($dist[$i]>$dist_min) {
      my $norm=1;
      if ( "$type" eq "bond" ) {
	$norm=$x[$i]*$x[$i];
      } elsif ( "$type" eq "angle" ) {
	$norm=sin($x[$i]);
      }
      $pot[$i]=-$kbT*log($dist[$i]/$norm);
    }
    else {
      $pot[$i]="nan";
      $flag[$i]="u";
    }
}

#find a valid point
my $valid_i=-1;
for (my $i=0;$i<=$#pot;$i++){
   if ($flag[$i] eq "i") {
     $valid_i=$i;
     last;
   }
}
die "All data points from file '$infile' are invalid after Boltzmann inversion, please check if your distribution is a valid dist.\n" if ($valid_i==-1);

#set point at beginning to invalid
my $first=undef;
for (my $i=$valid_i;$i>=0;$i--){
   if ($flag[$i] eq "u") {
     $pot[$i]=$pot[$i+1];
     $flag[$i]="o";
     $first=$i unless defined $first;
   }
}
$first=0 unless defined $first;

#set point at end to invalid
my $last=undef;
for (my $i=$valid_i;$i<=$#pot;$i++){
   if ($flag[$i] eq "u") {
     $pot[$i]=$pot[$i-1];
     $flag[$i]="o";
     $last=$i unless defined $last;
   }
}
$last=$#pot unless defined $last;
my $n=10;
my $valid=$last-$first+1;
die "Only $valid points are valid after Boltzmann inversion from file '$infile', please check if your distribution is a valid dist.\n" if ($valid < $n);

saveto_table($outfile,@x,@pot,@flag) || die "$progname: error at save table\n";
