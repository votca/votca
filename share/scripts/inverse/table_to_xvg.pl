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
my $type="non-bonded";
my $gmx_max=undef;

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

This script converts csg potential files to the xvg format.

Allowed options:
-v, --version         print version
-h, --help            show this help message
--type XXX            change the type of xvg table
                      Default: $type
--max MAX             Replace all pot value bigger MAX by MAX 


Possible types: non-bonded (=C12), bond, thermoforce, C12, C6

Examples:
* $progname --type bond table.in table_b0.xvg
EOF
    exit 0;
  }
  elsif ($ARGV[0] eq "--type"){
      shift(@ARGV);
      $type = shift(@ARGV);
  }
  elsif ($ARGV[0] eq "--max"){
      shift(@ARGV);
      $gmx_max = shift(@ARGV);
  }
  else{
    die "Unknow option '".$ARGV[0]."' !\n";
  }
}

die "2 parameters are nessary\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my @r;
my @pot;
my @flag;
(readin_table($infile,@r,@pot,@flag)) || die "$progname: error at readin_table\n";

if (defined($gmx_max)) {
  #gromacs does not like VERY big numbers
  for (my $i=0;$i<=$#r;$i++) {
    $pot[$i]=$gmx_max if $pot[$i]>$gmx_max;
    $pot[$i]=-$gmx_max if $pot[$i]<-$gmx_max;
  }
}

my @force;

#calc force
for (my $i=1;$i<$#r;$i++){
   $force[$i]=-($pot[$i+1]-$pot[$i-1])/($r[$i+1]-$r[$i-1]);
}
if ( "$type" eq "dihedral" ) {
  $force[0]=-($pot[1]-$pot[$#r-1])/($r[1]-$r[0]+$r[$#r]-$r[$#r-1]);
  $force[$#r]=$force[0];
} else {
  $force[0]=0;
  $force[$#r]=0.0;
}

open(OUTFILE,"> $outfile") or die "saveto_table: could not open $outfile\n";

#preserve comments
open(INFILE, "$infile");
while (<INFILE>){
	if($_ =~ /^[#@]/){print OUTFILE $_;}
}
close(INFILE);

my $fmt=undef;
my $begin=0;
my $end=undef;
if (( "$type" eq "non-bonded" ) or ("$type" eq "C12" )) {
  $fmt=sprintf("%%15.10e   %15.10e %15.10e   %15.10e %15.10e   %%15.10e %%15.10e\n",0,0,0,0);
}
elsif ( "$type" eq "C6" ){
  $fmt=sprintf("%%15.10e   %15.10e %15.10e   %%15.10e %%15.10e   %15.10e %15.10e\n",0,0,0,0);
}
elsif ( "$type" eq "bonded" ){
  $fmt="%15.10e   %15.10e %15.10e\n";
}
elsif ( "$type" eq "angle" ){
  $fmt="%15.10e   %15.10e %15.10e\n";
  $end=180;
}
elsif ( "$type" eq "dihedral" ){
  $fmt="%15.10e   %15.10e %15.10e\n";
  $begin=-180;
  $end=180;
}
elsif ( "$type" eq "thermforce" ){
  $fmt="%15.10e   %15.10e %15.10e\n";
}
else{
  die "$progname: Unsupported type of interatction: $type -> go and implement it\n";
}

die "$progname: table for type $type should begin with $begin, but I found $r[0]\n" if(abs($begin-$r[0]) > 1e-3);
die "$progname: table for type $type should end with $end, but I found $r[$#r]\n" if(($end) and (abs($end-$r[$#r]) > 1e-3));

for(my $i=0;$i<=$#r;$i++){
    printf(OUTFILE "$fmt",$r[$i],$pot[$i], $force[$i]);
}
close(OUTFILE) or die "Error at closing $outfile\n";

