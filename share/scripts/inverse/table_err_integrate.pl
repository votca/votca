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
	else
	{
		die "Unknow option '".$ARGV[0]."' !\n";
	}
}

#Print usage
die "no files given\n$usage\n" unless $#ARGV > 0;

use CsgFunctions;

sub readin_force_table($\@\@\@) {
  defined($_[3]) || die "readin_force_table: Missing argument\n";
  open(TAB,"$_[0]") || die "readin_force_table: could not open file $_[0]\n";
  my $line=0;
  while (<TAB>){
    $line++;
    # remove leading spacees for split
    $_ =~ s/^\s*//;
    next if /^[#@]/;
    next if /^\s*$/;
    my @parts=split(/\s+/);
    defined($parts[2]) || die "readin_force_table: Not enough columns in line $line in file $_[0]\n";
    #very trick array dereference (@) of pointer to an array $_[.] stored in an array $_
    push(@{$_[1]},$parts[0]);
    push(@{$_[2]},$parts[1]);
    push(@{$_[3]},$parts[2]);
  }
  close(TAB) || die "readin_force_table: could not close file $_[0]\n";
  return $line;
}

sub saveto_force_table($\@\@\@) {
  defined($_[3]) || die "saveto_force_table: Missing argument\n";
  open(OUTFILE,"> $_[0]") or die "saveto_force_table: could not open $_[0] \n";
  for(my $i=0;$i<=$#{$_[1]};$i++){
    print OUTFILE "${$_[1]}[$i] ${$_[2]}[$i] ${$_[3]}[$i]\n";
  }
  close(OUTFILE) or die "Error at closing $_[0]\n";
  return 1;
}

my $infile="$ARGV[0]";
my @r;
my @force;
my @force_errors;

readin_force_table($infile,@r,@force,@force_errors) || die "$progname: error at readin_force_table\n";

my $outfile="$ARGV[1]";
my @pot;
my @pot_errors;
my @ww;

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

saveto_force_table($outfile,@r,@pot,@pot_errors) || die "$progname: error at save table\n";
