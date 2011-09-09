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
my $epsilon=1e-5;
my $output=undef;
my $weak=undef;

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
    if($ARGV[0] eq "--output") {
        shift(@ARGV);
        $output = shift(@ARGV);
    }
    elsif($ARGV[0] eq "--epsilon") {
        shift(@ARGV);
        $epsilon = shift(@ARGV);
    }
    elsif($ARGV[0] eq "--weak") {
        shift(@ARGV);
	$weak="yes";
    }
    elsif (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help"))
	{
            print <<EOF;
$progname, version %version%
This script compares two tables

$usage

Allowed options:
    --output FILE     Output differnce to a file (used if --weak)
    --weak            Do not die, if y values or flags differ, just print sqrt difference.
    --error  ERR      Relative error
                      Default: $epsilon 
-h, --help            Show this help message
EOF
  exit 0;
    }
}

die "2 parameters are nessary\n" if ($#ARGV<1);

use CsgFunctions;

my $file1="$ARGV[0]";
my $file2="$ARGV[1]";

my @r1;
my @pot1;
my @flag1;
(readin_table($file1,@r1,@pot1,@flag1)) || die "$progname: error at readin_table\n";

my @r2;
my @pot2;
my @flag2;
(readin_table($file2,@r2,@pot2,@flag2)) || die "$progname: error at readin_table\n";

$#r1 == $#r2 || die "$progname: error, tables have different length";

sub difference($$$) {
  defined($_[2]) || die "difference: Missing argument\n";
  my $x=$_[0];
  my $y=$_[1];
  $x=0 if ($x =~ /nan/i);
  $y=0 if ($y =~ /nan/i);
  my $error=$_[2];
  return 0 if $x == $y; #fetch the case that both are zero
  return abs($x-$y) if (abs($x-$y) < $error);
  return abs($x-$y)/abs($x) if (abs($x) > abs($y));
  return abs($x-$y)/abs($y);
}

my $sum=0;
for (my $i=0;$i<=$#r1; $i++) {
  abs($r1[$i] - $r2[$i]) < $epsilon || die "$progname: first column different at position $i\n";
  #check relative error!
  my $diff=&difference($pot1[$i],$pot2[$i],$epsilon);
  if ($weak) {
    $sum+=$diff;
  } else {
    die "$progname: second column different at position $i\n" if ($diff > $epsilon);
    $flag1[$i] eq $flag2[$i] || die "$progname: flag different at position $i\n";
  }
}
if ($output) {
  open(FILE,"> $output") or die "Could not open $output for write\n";
  print FILE "$sum";
  close(FILE) or die "Error at closing $output\n";
} else {
  print "$sum";
}

