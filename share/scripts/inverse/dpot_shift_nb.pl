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
use strict;

( my $progname = $0 ) =~ s#^.*/##;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
$progname, version %version%
This script set the beginning of the dpot to the first valid value and shift the whole potential
so that dpot(r_max)=0.

Usage: $progname infile outfile

NEEDS:

USES: readin_table saveto_table
EOF
  exit 0;
}

die "2 parameters are nessary, <infile> <outfile>\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

# read in the current dpot
my @r;
my @dpot;
my @flag;
(readin_table($infile,@r,@dpot,@flag)) || die "$progname: error at readin_table\n";

# find first u
my $i_first;
for($i_first=0; ($i_first<=$#r) && ($flag[$i_first] =~ /[u]/); $i_first++) {}

# shift beginning
for(my $i=0; $i<$i_first; $i++) {
    $dpot[$i] = $dpot[$i_first];
    $flag[$i]="o"
}

# bring end to zero
for(my $i=0; $i<=$#r; $i++) {
    $dpot[$i] -= $dpot[$#r];
}

# save to file
saveto_table($outfile,@r,@dpot,@flag) || die "$progname: error at save table\n";
