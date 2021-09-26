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

$_=$0;
s#^.*/##;
my $progname=$_;
my $usage="Usage: $progname [OPTIONS] <in> <out>";

my $with_errors="no";
my $with_entropie="no";
my $kbT=undef;
my $from="right";
my $spherical="no";

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
This script adjusts a distribution in such a way that value smaller 0 will be replaces with 0.

$usage

Allowed options:
-h, --help            Show this help message

Examples:
* $progname CG-CG.dist.tmp CG-CG.dist.new
END
		exit;
	}
	else
	{
		die "Unknown option '".$ARGV[0]."' !\n";
	}
}

die "2 parameters are necessary\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my @r;
my @dist;
my @flag;
my $comments="";

(readin_table($infile,@r,@dist,@flag,$comments)) || die "$progname: error at readin_table\n";
my $outfile="$ARGV[1]";
for (my $i=0;$i<=$#r;$i++){
  $dist[$i]=0 if ($dist[$i]<0);
}
(saveto_table($outfile,@r,@dist,@flag,$comments)) || die "$progname: error at save table\n";
