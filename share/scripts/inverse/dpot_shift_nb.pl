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

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
$progname, version %version%
This script shifts the whole potential to the last value, like it is normally done for non-bonded potentials.

Usage: $progname infile outfile
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
my $comments;
(readin_table($infile,@r,@dpot,@flag,$comments)) || die "$progname: error at readin_table\n";

# bring end to zero
for(my $i=0; $i<=$#r; $i++) {
    $dpot[$i] -= $dpot[$#r];
}

# save to file
saveto_table($outfile,@r,@dpot,@flag,$comments) || die "$progname: error at save table\n";
