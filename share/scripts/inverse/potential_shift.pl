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
my $usage="Usage: $progname [OPTIONS] <in> <out>";
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
  if (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help")){
    print <<EOF;
$progname, version %version%

This script shifts the whole potential by minimum (bonded potentials) or last value (non-bonded potentials).

$usage

Allowed options:
-h, --help            show this help message
--type XXX            change the type of potential
                      Default: $type

Possible types: non-bonded, bond, angle, dihedral, bonded

Examples:
* $progname --type bond table.in table.out
EOF
    exit 0;
  }
  elsif ($ARGV[0] eq "--type"){
      shift(@ARGV);
      $type = shift(@ARGV);
  }
  else{
    die "Unknown option '".$ARGV[0]."' !\n";
  }
}

die "2 parameters are necessary, <infile> <outfile>\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

# read in the current dpot
my @r;
my @dpot;
my @flag;
my $comments;
(readin_table($infile,@r,@dpot,@flag,$comments)) || die "$progname: error at readin_table\n";

my $zero=undef;
if (( "$type" eq "non-bonded" ) or ("$type" eq "thermforce" )) {
  $zero=$dpot[$#r];
} 
elsif (( "$type" eq "bond" ) or ("$type" eq "dihedral") or ("$type" eq "angle") or ("$type" eq "bonded")) {
  for(my $i=0; $i<=$#r; $i++) {
    $zero=$dpot[$i] if (($flag[$i] =~ /[i]/) and not defined($zero));
    $zero=$dpot[$i] if (($flag[$i] =~ /[i]/) and ($dpot[$i]<$zero));
  }
  die "No valid value found in $infile" unless defined($zero);
}
else{
  die "$progname: Unsupported type of interatction: $type -> go and implement it\n";
}

# shift potential by $zero
for(my $i=0; $i<=$#r; $i++) {
    $dpot[$i] -= $zero;
}

# save to file
saveto_table($outfile,@r,@dpot,@flag,$comments) || die "$progname: error at save table\n";
