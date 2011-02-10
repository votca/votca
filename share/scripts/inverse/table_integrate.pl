#! /usr/bin/perl -w
#
# Copyright 2009,2010 The VOTCA Development Team (http://www.votca.org)
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
my $usage="Usage: $progname [OPTIONS] <in> <out>";

my $with_errors="no";
my $with_entropie="no";
my $kbT=undef;

# read program arguments

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
This script calculates the integral of a table. Please note the force is the NEGATIVE integral of the potential (use 'table linearop' and multiply the table with -1)

$usage

Allowed options:
    --with-errors     calculate error
    --with-S          Add entropic contribution to force
                      2*k_B T/r
    --kbT   NUMBER    use NUMBER as k_B*t for the entropic part
-h, --help            Show this help message

Examples:
* $progname --with-S --kbT 2.49435 tmp.force tmp.dpot


USES: readin_table saveto_table
END
		exit;
	}
	elsif ($ARGV[0] eq "--with-errors"){
          shift(@ARGV);
	  $with_errors="yes";
	}
	elsif ($ARGV[0] eq "--with-S"){
          shift(@ARGV);
	  $with_entropie="yes";
	}
	elsif ($ARGV[0] eq "--kbT"){
          shift(@ARGV);
	  $kbT=shift(@ARGV);
	}
	else
	{
		die "Unknow option '".$ARGV[0]."' !\n";
	}
}

#Print usage
die "no files given\n$usage\n" unless $#ARGV > 0;

use CsgFunctions;

my $infile="$ARGV[0]";
my @r;
my @force;
my @flag;
my @force_errors;

if ("$with_errors" eq "yes") {
  (readin_table_err($infile,@r,@force,@force_errors,@flag)) || die "$progname: error at readin_table\n";
} else {
  (readin_table($infile,@r,@force,@flag)) || die "$progname: error at readin_table\n";
}

if ("$with_entropie" eq "yes"){
  die "$progname: kbT not defined specify it with --kbT option\n" unless defined($kbT);
  for (my $i=$#r;$i>=0;$i--){
    if ($r[$i]>0) {
      $force[$i] += 2*$kbT/$r[$i];
    }
  }
}

my $outfile="$ARGV[1]";
my @pot;
my @pot_errors;
my @ww;

#calc pot with trapez rule
#int_j= sum_i^j (r_i+1 - r_i)*(f_i+f_i+1)/2
#int_j+1= int_j + (r_i+1 - r_i)*(f_i+f_i+1)/2
#int_j= int_j+1 - (r_i+1 - r_i)*(f_i+f_i+1)/2
#begin from end to make pot(max)=0
$pot[$#r]=0;
$ww[$#r]=0;
for (my $i=$#r-1;$i>=0;$i--){
  #hh = delta x /2
  my $hh=0.5*($r[$i+1]-$r[$i]);
  $pot[$i]=$pot[$i+1] - $hh*($force[$i+1]+$force[$i]);
  $ww[$i]+= $hh;
  $ww[$i+1]+= $hh;
}
#ww contains w_i=(r_i+1-r_i-1)/2


if ("$with_errors" eq "yes") {
  #all error are independent(we assume that)
  #resort sum (only one force per summand)
  # U_j= sum_i ^j = sum_i^j f_i(r_i+1 - r_i-1)/2 + randterm
  # o^2(U_j)=sum_i o^2(f_i)*(r_i+1 - r_i-1)/2 + o^2(randterm)
  my $var_int = ($ww[$#r]*$force_errors[$#r])**2;
  $pot_errors[$#r]=sqrt($var_int);
  for(my $i=$#r-1; $i>=0;$i--) {
    my $hh = 0.5*($r[$i+1] - $r[$i]);
    $pot_errors[$i] = sqrt($var_int + ($hh*$force_errors[$i])**2);
    $var_int += ($ww[$i]*$force_errors[$i])**2;
  }
}

if ("$with_errors" eq "yes") {
  saveto_table_err($outfile,@r,@pot,@pot_errors,@flag) || die "$progname: error at save table\n";
}else {
  saveto_table($outfile,@r,@pot,@flag) || die "$progname: error at save table\n";
}
