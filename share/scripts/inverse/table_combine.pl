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
my $usage="Usage: $progname [OPTIONS] <in> <in2> <out>";
my $epsilon=1e-5;
my $op=undef;
my $noflags=undef;
my $dosum=undef;
my $die=undef;
my $scale=1.0;

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
    elsif($ARGV[0] eq "--error") {
        shift(@ARGV);
        $epsilon = shift(@ARGV);
    }
    elsif($ARGV[0] eq "--no-flags") {
        shift(@ARGV);
	$noflags="yes";
    }
    elsif($ARGV[0] eq "--op") {
        shift(@ARGV);
	$op = shift(@ARGV);
    }
    elsif($ARGV[0] eq "--scale") {
        shift(@ARGV);
	$scale = shift(@ARGV);
    }
    elsif($ARGV[0] eq "--die") {
        shift(@ARGV);
	$die = "yes";
    }
    elsif($ARGV[0] eq "--sum") {
        shift(@ARGV);
	$dosum = "yes";
    }
    elsif (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help"))
	{
            print <<EOF;
$progname, version %version%
This script combines two tables with a certain operation

$usage

Allowed options:
    --error  ERR      Relative error
                      Default: $epsilon
    --op OP           Operation to perform
                      Possible: =,+,-,,/,d,x
		      d = |y1-y2|, x=* (to avoid shell trouble)
    --sum             Output the sum instead of a new table
    --die             Die if op '=' fails
    --no-flags        Do not check for the flags
    --scale XXX       Scale output/sum with this number
                      Default $scale
-h, --help            Show this help message
EOF
  exit 0;
    }
    else
    {
      die "Unknown option '".$ARGV[0]."' !\n";
    }
}

die "$progname: Missing --op options" unless ($op);

use CsgFunctions;

if ($die || $dosum) {
  die "2 parameters are necessary\n" if ($#ARGV<1);
} else {
  die "3 parameters are necessary\n" if ($#ARGV<2);
}
my $file1="$ARGV[0]";
my $file2="$ARGV[1]";

my @r1;
my @pot1;
my @flag1;
my $comments1;
(readin_table($file1,@r1,@pot1,@flag1,$comments1)) || die "$progname: error at readin_table\n";

my @r2;
my @pot2;
my @flag2;
my $comments2;
(readin_table($file2,@r2,@pot2,@flag2,$comments2)) || die "$progname: error at readin_table\n";

$#r1 == $#r2 || die "$progname: error, tables have different length";

sub difference_relative($$$) {
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

sub operation($$$) {
  defined($_[2]) || die "operation: Missing argument\n";
  my $x=$_[0];
  my $op="$_[1]";
  my $y=$_[2];
  use Switch;
  switch($op) {
    case /\+|-|\*|\/|x/ {
      $op="*" if ($op eq "x");
      my $val = eval "$x $op $y";
      die "operation: Could not calculate '$x $op $y'\n" if $@;
      return $val;
    }
    case "=" {
      my $diff=&difference_relative($x,$y,$epsilon);
      return 1 if ($diff > $epsilon);
      return 0;
    }
    case "d" {
      return abs($x-$y);
    } 
    else {
      die "operation: Unknown operation $op\n";
    }
  }
}

my $sum=0;
my @table;
for (my $i=0;$i<=$#r1; $i++) {
  abs($r1[$i] - $r2[$i]) < $epsilon || die "$progname: first column different at position $i\n";
  my $value=&operation($pot1[$i],$op,$pot2[$i]);
  if (($die)&&($op eq "=")&&($value == 1)) {
    die "progname: second column different at position $i\n";
  }
  $value*=$scale;
  $sum+=$value;
  $table[$i]=$value;
  unless ($noflags){
    $flag1[$i] eq $flag2[$i] || die "$progname: flag different at position $i\n";
  }
}

if ($die) {
  #notthing
} elsif ($dosum) {
  print "$sum\n";
} else {
  my $comments="# $progname: combining $file1 $op $file2 into $ARGV[2]\n";
  $comments.="#Comments from $file1\n$comments1" if (defined($comments1));
  $comments.="#Comments from $file2\n$comments2" if (defined($comments2));
  saveto_table($ARGV[2],@r1,@table,@flag1,$comments) || die "$progname: error at save table\n";
}
