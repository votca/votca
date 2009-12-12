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
#
# extrapolation methods:
#
# linear:
#   y = ax + b;  b = - m*xp + yp; a = m
# sasha_shit:
#   y = a*(x-b)^2; b = (x0 - 2y0/m); a = m^2/(4*y0)
# exponential:
#   y = a*exp(b*x); a = y0*exp(-m*x0/y0); b = m/y0;
# quadratic:
#   y = curv*(x+a)^2 + b; a = m/(2*curv) - x0; b = y0 - m^2/(4*curv)
#

sub extrapolate_linear($$$$) {
  my $x0 = $_[0];
  my $y0 = $_[1]; 
  my $m =  $_[2];
  my $x =  $_[3];

  return $m*($x - $x0) + $y0;
}

sub sasha_shit($$$$) {
  my $x0 = $_[0];
  my $y0 = $_[1]; 
  my $m =  $_[2];
  my $x =  $_[3];

  my $a = ($m**2)/(4*$y0);
  my $b = $x0 - 2*$y0/$m;
  #my $b = $x0 + 2*$y0/$m;
  #my $a = $m/(2*($x0-$b));

  return $a*($x-$b)**2;
}
sub extrapolate_quad($$$$) {
  my $x0 = $_[0];
  my $y0 = $_[1]; 
  my $m =  $_[2];
  my $x =  $_[3];

  # $curv is a global variable
  my $a = 0.5*$m/$curv - $x0;
  my $b = $y0 - 0.25*$m*$m/$curv;

  return $curv*($x + $a)**2 + $b
}
sub extrapolate_exp($$$$) {
  my $x0 = $_[0];
  my $y0 = $_[1]; 
  my $m =  $_[2];
  my $x =  $_[3];

  my $a = $y0*exp(-$m*$x0 / $y0);
  my $b = $m/$y0;

  return $a*exp($b*$x);
}

use strict;

$_=$0;
s#^.*/##;
my $progname=$_;
my $usage="Usage: $progname [OPTIONS] <in> <out>";

my $avgpoints = 3;
my $function="quadratic";
my $region = "leftright";
our $curv = 10000.0; # curvature for quadratic extrapolation

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
    if($ARGV[0] eq "--avgpoints") {
        $avgpoints = $ARGV[1];
        shift(@ARGV);
        shift(@ARGV);
    }
    elsif($ARGV[0] eq "--function") {
        $function = $ARGV[1];
        shift(@ARGV);
        shift(@ARGV);
    }
    elsif($ARGV[0] eq "--region") {
        $region = $ARGV[1];
        shift(@ARGV);
        shift(@ARGV);
    }
    elsif($ARGV[0] eq "--curvature") {
        $curv = $ARGV[1];
        shift(@ARGV);
        shift(@ARGV);
    }
    elsif (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help"))
	{
		print <<END;
$progname, version %version%
This script extrapolates a table

$usage

Allowed options:
--avgpoints           average over so many points to extrapolate: default is 3
--function            linear, quadratic or exponential, sasha: default is quadratic
--region              left, right, or leftright: default is leftright
--curvature           curvature of the quadratic function: default is 10000
                      makes sense only for quadratic extrapolation, ignored for other cases
-h, --help            Show this help message

NEEDS:

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

my $infile="$ARGV[0]";
my @r;
my @val;
my @flag;
(readin_table($infile,@r,@val,@flag)) || die "$progname: error at readin_table\n";

my $outfile="$ARGV[1]";

#==============
my ($do_left, $do_right);

# parse $region: decide where to extrapolate
if ($region eq "left") {
  $do_left = 1;
  $do_right = 0;
}
elsif ($region eq "right") {
  $do_left = 0;
  $do_right = 1;
}
elsif ($region eq "leftright") {
  $do_left = 1;
  $do_right = 1;
}
else {
  die "$progname: Unknown region: $region !\n";
}

my $extrap_method;

# parse $function: decide which method to use
if ($function eq "linear") {
   $extrap_method = \&extrapolate_linear;
}
elsif ($function eq "quadratic") {
   $extrap_method = \&extrapolate_quad;
}
elsif ($function eq "exponential") {
   $extrap_method = \&extrapolate_exp;
}
elsif ($function eq "sasha") {
   $extrap_method = \&sasha_shit;
}
else {
  die "$progname: Unknown extrapolation function: $function !\n";
}

# do extrapolation: left
if ($do_left) {
  # find beginning
  my $first;
  for ($first=1;$first<=$#r;$first++) {
     last if($flag[$first] eq "i");
  }

  # grad  of beginning
  my $grad_beg;
  $grad_beg = ($val[$first + $avgpoints] - $val[$first])/($r[$first + $avgpoints] - $r[$first]);

  # now extrapolate beginning
  for(my $i=$first-1; $i >= 0; $i--) {
      $val[$i] = &{$extrap_method}($r[$first], $val[$first], $grad_beg, $r[$i]);
  }
}

# do extrapolation: right
if ($do_right) {
  # find end
  my $last;
  for ($last=$#r;$last>0;$last--) {
     last if($flag[$last] eq "i");
  }

  # grad  of end
  my $grad_end;
  $grad_end = ($val[$last] - $val[$last - $avgpoints])/($r[$last] - $r[$last-$avgpoints]);

  # now extrapolate ends
  for(my $i=$last+1; $i <= $#r; $i++) {
      $val[$i] = &{$extrap_method}($r[$last], $val[$last], $grad_end, $r[$i]);
  } 
}
#==============


saveto_table($outfile,@r,@val,@flag) || die "$progname: error at save table\n";
