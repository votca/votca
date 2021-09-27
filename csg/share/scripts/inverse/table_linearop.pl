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

$_=$0;
s#^.*/##;
my $progname=$_;
my $usage="Usage: $progname [OPTIONS] <in> <out> <a> <b>";

#Defaults
my $withflag=undef;
my $with_errors="no";
my $col="y";

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
This script performs a linear operation on the y values:
''\$y_{new} = a*y_{old} + b\$''

$usage

Allowed options:
-h, --help            Show this help message
    --withflag  FL    only change entries with specific flag in src
    --with-errors     also read and calculate errors
    --on-x            work on x values instead of y values

Examples:
* $progname tmp.dpot.cur tmp.dpot.new 1.0 0.0
END
		exit;
	}
    elsif ($ARGV[0] eq "--withflag"){
        shift(@ARGV);
        die "nothing given for --withflag" unless $#ARGV > -1;
        $withflag = shift(@ARGV);
    }
    elsif ($ARGV[0] eq "--with-errors"){
          shift(@ARGV);
	  $with_errors="yes";
    }	
    elsif ($ARGV[0] eq "--on-x"){
          shift(@ARGV);
	  $col="x";
    }	
    else
	{
		die "Unknown option '".$ARGV[0]."' !\n";
	}
}

#Print usage
die "missing parameters\n$usage\n" unless $#ARGV >= 3;

my $a = $ARGV[2];
my $b = $ARGV[3];

use CsgFunctions;

my $file="$ARGV[0]";
my $outfile="$ARGV[1]";

print "$progname: $file to $outfile with $col' = $a*$col + $b \n";

my @r;
my @val;
my @flag;
my @errors;
my $comments="";
if ("$with_errors" eq "yes") {
  (readin_table_err($file,@r,@val,@errors,@flag,$comments)) || die "$progname: error at readin_table\n";
} else {
  (readin_table($file,@r,@val,@flag,$comments)) || die "$progname: error at readin_table\n";
}

for(my $i=0; $i<=$#r; $i++) {
  # skip if flag does not match
  if(($withflag) and ($flag[$i] !~ m/[$withflag]/)) {
      next;
  }
  if ("$col" eq "x") {
    $r[$i]=$a*$r[$i]+$b;
  } else {
    $val[$i] = $a*$val[$i] + $b;
    if ("$with_errors" eq "yes") {
      $errors[$i] = $a*$errors[$i];
    }
  }
}

$comments.="# $progname: $file -> $outfile $col' = $a*$col + $b\n";
if ("$with_errors" eq "yes") {
  saveto_table_err($outfile,@r,@val,@errors,@flag,$comments) || die "$progname: error at save table\n";
}else {
  saveto_table($outfile,@r,@val,@flag,$comments) || die "$progname: error at save table\n";
}
