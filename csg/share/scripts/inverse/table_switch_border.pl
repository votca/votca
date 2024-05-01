#! /usr/bin/env -S perl -w
#
# Copyright 2016 The VOTCA Development Team (http://www.votca.org)
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
This script applies a switching function to the end of the table to switch it smoothly to zero by y = y*cos( pi*(x-x_switch)/(2*(x_end-x_switch)) )

Usage: $progname infile outfile <x_switch>
EOF
  exit 0;
}

die "3 parameters are necessary\n" if ($#ARGV<2);

use CsgFunctions;

my $infile="$ARGV[0]";
my @r_cur;
my @pot_cur;
my @flag_cur;
(readin_table($infile,@r_cur,@pot_cur,@flag_cur)) || die "$progname: error at readin_table\n";

my $outfile="$ARGV[1]";
my @pot;
my $a = $ARGV[2];

# TODO: think about addition rules
# now I did it like that to always maintain interval of interest in all potentials

# find end
my $last;
for ($last=$#r_cur;$last>0;$last--) {
   last if($flag_cur[$last] eq "i");
}

use constant { 
   PI => 4 * atan2(1,1) 
};

for (my $i=0;$i<=$#r_cur;$i++){  
  $pot[$i]=$pot_cur[$i];
  if($flag_cur[$i] eq "i") {
    if($r_cur[$i]>$a) {
        $pot[$i] = $pot_cur[$i] * cos(PI*($r_cur[$i]-$a)/(2.0*($r_cur[$last]-$a)))    
    }
  }
}


saveto_table($outfile,@r_cur,@pot,@flag_cur) || die "$progname: error at save table\n";

