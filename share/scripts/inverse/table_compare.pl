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

if (defined($ARGV[0]) && ( ("$ARGV[0]" eq "-h" ) || ("$ARGV[0]" eq "--help") )){
  print <<EOF;
$progname, version %version%
This script compares two tables

Usage: $progname infile1 infile2
EOF
  exit 0;
}

die "2 parameters are nessary\n" if ($#ARGV<1);

use CsgFunctions;

my $epsilon=1e-5;

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

for (my $i=0;$i<=$#r1; $i++) {
  abs($r1[$i] - $r2[$i]) < $epsilon || die "$progname: first column different at position $i\n";
  #check relative error!
  abs($pot1[$i] - $pot2[$i])/(($pot1[$i] == 0.0) ? 1.0 : $pot1[$i]) < $epsilon || die "$progname: second column different at position $i\n";
  $flag1[$i] eq $flag2[$i] || die "$progname: flag different at position $i\n";
}
