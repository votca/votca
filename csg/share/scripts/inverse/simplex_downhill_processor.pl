#! /usr/bin/env -S perl -w
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
Changes a simplex state according to the current state using the Nelder–Mead method or downhill simplex algorithm.

Usage: $progname current_state new_state
EOF
  exit 0;
}

die "2 parameter are necessary\n" if ($#ARGV<1);

use CsgSimplex;

my $alpha=1; #Reflection constant
my $gamma=2; #Expansion constant
my $rho=0.5; #Contraction constant
my $sigma=0.5; #Reduction constant

my @simplex_table;
my $state;
my $comments;
my $parameter_names="not in state file";
(readin_simplex_state($ARGV[0],$state,@simplex_table,$comments,$parameter_names)) || die "$progname: error at readin_simplex_table\n";
print "We are in state $state with parameters ($parameter_names):\n";
for (my $i=0;$i<=$#simplex_table;$i++){
  print "@{$simplex_table[$i]}\n";
}

my $highest=get_convergence_value(@simplex_table,"highest");
my $second_highest=get_convergence_value(@simplex_table,"second");
my $lowest=get_convergence_value(@simplex_table,"lowest");
my $try=get_convergence_value(@simplex_table,"try");
print "values: $highest (highest), $second_highest (2nd highest), $lowest (lowest), $try (try)\n";
my $next_state;
if ($state eq "Initialization") {
  replace_parameter_flag(@simplex_table,"try","complete");
  $next_state="Reflection";
} elsif ($state eq "Reflection") {
  if ($try < $lowest) {
    $next_state="Expansion";
    #we need to remember this value
    replace_parameter_flag(@simplex_table,"try","tryold");
  } elsif ( $try > $second_highest ) {
    $next_state="Contraction";
    if ($try > $highest ) {
      remove_parameter_set(@simplex_table,"try");
    } else {
      replace_parameter_flag(@simplex_table,"try","complete");
      remove_parameter_set(@simplex_table,"highest");
    }
  } else { #$try is between $lowest and $second_highest
    $next_state="Reflection";
    replace_parameter_flag(@simplex_table,"try","complete");
    remove_parameter_set(@simplex_table,"highest");
  }
} elsif ($state eq "Expansion") {
  $next_state="Reflection";
  my $tryold=get_convergence_value(@simplex_table,"tryold");
  if ($try < $tryold) { #tryold is the reflection point from before
    remove_parameter_set(@simplex_table,"tryold");
    replace_parameter_flag(@simplex_table,"try","complete");
    remove_parameter_set(@simplex_table,"highest");
  } else {
    remove_parameter_set(@simplex_table,"try");
    replace_parameter_flag(@simplex_table,"tryold","complete");
    remove_parameter_set(@simplex_table,"highest");
  }
} elsif ($state eq "Contraction") {
  if ($try < $highest) {
    replace_parameter_flag(@simplex_table,"try","complete");
    remove_parameter_set(@simplex_table,"highest");
    $next_state="Reflection";
  } else {
    $next_state="Reduction";
    remove_parameter_set(@simplex_table,"try");
  }
} elsif ($state eq "Reduction") {
  replace_parameter_flag(@simplex_table,"try","complete");
  $next_state="Reflection";
}else {
  die "$progname: Unknown state '$state'\n";
}

if ($next_state eq "Reflection") {
  sort_simplex_table(@simplex_table);
  my @center_parameter=calc_parameter_center(@simplex_table);
  my @highest_parameter=@{$simplex_table[$#simplex_table]};
  my @try_paramter=linop_parameter(@center_parameter,$alpha,@center_parameter,@highest_parameter);
  push(@simplex_table,\@try_paramter);
} elsif ($next_state eq "Expansion") {
  my @tryold_parameter=remove_parameter_set(@simplex_table,"tryold"); #this should not go into the center
  sort_simplex_table(@simplex_table);
  my @center_parameter=calc_parameter_center(@simplex_table);
  my @highest_parameter=@{$simplex_table[$#simplex_table]};
  my @try_paramter=linop_parameter(@center_parameter,$gamma,@center_parameter,@highest_parameter);
  push(@simplex_table,\@try_paramter,\@tryold_parameter);
} elsif ($next_state eq "Contraction") {
  sort_simplex_table(@simplex_table);
  my @center_parameter=calc_parameter_center(@simplex_table);
  my @highest_parameter=@{$simplex_table[$#simplex_table]};
  my @try_paramter=linop_parameter(@highest_parameter,$rho,@center_parameter,@highest_parameter);
  push(@simplex_table,\@try_paramter);
} elsif ($next_state eq "Reduction") {
  sort_simplex_table(@simplex_table);
  my @lowest_parameter=@{$simplex_table[0]};
  for (my $i=1; $i<=$#simplex_table;$i++) {
    my @try_paramter=linop_parameter(@lowest_parameter,$sigma,@{$simplex_table[$i]},@lowest_parameter);
    $simplex_table[$i]=\@try_paramter;
  }
}else {
  die "$progname: Unknown state '$next_state'\n";
}

print "Preparing $next_state with parameters ($parameter_names):\n";
for (my $i=0;$i<=$#simplex_table;$i++){
  print "@{$simplex_table[$i]}\n";
}
(saveto_simplex_state($ARGV[1],$next_state,@simplex_table,$comments)) || die "$progname: error at readin_simplex_table\n";
