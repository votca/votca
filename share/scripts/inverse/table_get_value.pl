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
my $usage="Usage: $progname [OPTIONS] X infile";

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
This script print the y value of x, which is closest to X.

$usage

Allowed options:
-h, --help            Show this help message

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

my $X="$ARGV[0]";
my $infile="$ARGV[1]";

my @x;
my @y;
my @flag;
(readin_table($infile,@x,@y,@flag)) || die "$progname: error at readin_table\n";

my $value=$y[0];
for(my $i=1; $i<=$#x; $i++) {
  if($x[$i]<$X) {
    $value=$y[$i];
  } else {
    $value=$y[$i] unless (($x[$i]-$X)>($X-$x[$i-1]));
    print "$value\n";
    exit 0;
  }
}
die "$progname: value $X not found\n";
