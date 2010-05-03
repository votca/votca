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

use strict;

( my $progname = $0 ) =~ s#^.*/##;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
$progname, %version%
This script calculates the potential from parameters
given in file simplex.new
In addition it does some magic tricks:
- does not crash when calc r(0)
- does not allow values of pot>1e10
- shifts the potential, so that it is zero at the cutoff

Usage: $progname infile outfile

USES: readin_table readin_simplex_table calc_func saveto_table

NEEDS: -

EOF
  exit 0;
}

die "3 parameters are nessary\n" if ($#ARGV<2);

use CsgFunctions;
use SimplexFunctions;

my $outfile="$ARGV[0]";
my $simplex_table="$ARGV[1]";
my $param_N="$ARGV[2]";
my $p_line_nr="$ARGV[3]";
my %hash=();

my @ftar;
my @params;
my @flag;

open(TAB,"$simplex_table") || die "could not open file $simplex_table\n";
my $line=0;
while (<TAB>) {
  $line++;
  $_ =~ s/^\s*//;
  next if /^[#@]/;
  next if /^\s*$/;
  my @values=split(/\s+/);
  defined($values[1]) || die "readin_table: Not enough columns in line $line in file $simplex_table\n";
  foreach (1..$param_N) {
    push @{$hash{"p_$_"}}, $values[$_];
  }
}
close(TAB) || die "could not close file $simplex_table\n";

@ftar=@{$hash{p_0}};
@flag=@{$hash{\"p_$param_N\"}};

# Mark current parameter set as active
$flag[$p_line_nr]="active";

saveto_simplex_table("simplex.new",@ftar,@sig,@eps,@flag_simplex) || die "$progname: error at saveto_simplex_table\n";