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

This script converts csg potential files to the xvg format.

In addition, it does some magic tricks:
- bigger value will be set to pot_max (see xml)
- all values after the cutoff are set to the same value as the cufoff

For non-bonded potential:
- shift the potential, so that it is zero at the cutoff
- put the values in the C12 column

Allowed options:
-v, --version         print version
-h, --help            show this help message
--type XXX            change the type of xvg table (non-bonded, bonded, ...)
                      Default: $type

Examples:
* $progname --type bonded table.in table_b0.xvg


USES: csg_get_property saveto_table readin_table
EOF
    exit 0;
  }
  elsif ($ARGV[0] eq "--type"){
      shift(@ARGV);
      $type = shift(@ARGV);
  }
  else{
    die "Unknow option '".$ARGV[0]."' !\n";
  }
}

die "2 parameters are nessary\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my $gromacs_max=csg_get_property("cg.inverse.gromacs.pot_max");
my $table_end=csg_get_property("cg.inverse.gromacs.table_end");
my $table_bins=csg_get_property("cg.inverse.gromacs.table_bins");

my @r;
my @pot;
my @flag;
(readin_table($infile,@r,@pot,@flag)) || die "$progname: error at readin_table\n";

#gromacs does not like VERY big numbers
for (my $i=0;$i<=$#r;$i++) {
  $pot[$i]=$gromacs_max if $pot[$i]>$gromacs_max;
  $pot[$i]=-$gromacs_max if $pot[$i]<-$gromacs_max;
}

#cutoff is last point
my $i_cut=$#r;

if ( "$type" eq "non-bonded" ) {
  #shift potential so that it is zero at cutoff
  for (my $i=0;$i<=$i_cut;$i++){
    $pot[$i]-=$pot[$i_cut];
  }
}

# set end of the potential to zero
for (my $i=$i_cut+1;$i<=$table_end/$table_bins;$i++) {
  $pot[$i]=$pot[$i_cut];
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

#preserve comments
open(INFILE, "$infile");
while (<INFILE>){
	if($_ =~ /^[#@]/){print OUTFILE $_;}
}
close(INFILE);

my $fmt=undef;
if ( "$type" eq "non-bonded" ){
  $fmt=sprintf("%%15.10e   %15.10e %15.10e   %15.10e %15.10e   %%15.10e %%15.10e\n",0,0,0,0);
}
elsif ( "$type" eq "bonded" ){
  $fmt="%15.10e   %15.10e %15.10e\n";
}
else{
  die "$progname: Unsupported type of interatction: $type -> go and implement it\n";
}

for(my $i=0;$i<=$#r;$i++){
    printf(OUTFILE "$fmt",$r[$i],$pot[$i], $force[$i]);
}
close(OUTFILE) or die "Error at closing $outfile\n";

