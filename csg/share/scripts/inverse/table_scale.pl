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
my $usage="Usage: $progname [OPTIONS] infile outfile prefactor1 prefactor2";

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
This script applies a prefactor to infile. The prefactor is is interpolated lines between the prefactor1 and prefactor2.
$usage
Allowed options:
-h, --help            Show this help message
END
		exit;
	}
	else
	{
		die "Unknown option '".$ARGV[0]."' !\n";
	}
}

#Print usage
die "no files given\n$usage\n" unless $#ARGV >=3 ;

use CsgFunctions;

my $infile="$ARGV[0]";
my @r;
my @val;
my @flag;
my $comments;
(readin_table($infile,@r,@val,@flag,$comments)) || die "$progname: error at readin_table\n";

my $outfile="$ARGV[1]";
my @out;

my $prefactor="$ARGV[2]";
my $prefactor2 = "$ARGV[3]";

for (my $i=0;$i<=$#r;$i++){
# do a linear interpoltation between the prefactors
   $out[$i]=$i/$#r*$val[$i]*$prefactor2+(1-$i/$#r)*$val[$i]*$prefactor;
}

saveto_table($outfile,@r,@out,@flag,$comments) || die "$progname: error at save table\n";
