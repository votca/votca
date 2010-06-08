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
$progname, version %version%
This script:
- Reads in current simplex table
- Reads in current state file
- Decides which transformation to perfom and saves it in new state file
- Adds new parameter set obtained from this transformation to simplex table

Usage: $progname infile outfile param_N

NEEDS: name ftol

USES: csg_get_property readin_simplex_table saveto_simplex_table 
EOF
  exit 0;
}

die "2 parameters are nessary\n" if ($#ARGV<1);

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
my $param_N="$ARGV[2]";

my $ndim=$param_N+1;
my $mdim=$ndim+1;

use CsgFunctions;
use SimplexFunctions;
use Switch;

my $name=csg_get_property("cg.non-bonded.name");
my $ftol=csg_get_property("cg.inverse.simplex.ftol");

# Read in simplex table and assign to arrays
my (%hash)=readin_simplex_table($infile,$ndim) or die "$progname: error at readin_simplex_table\n";
my @ftar=@{$hash{p_0}};
my @flag=@{$hash{"p_$ndim"}};

my @p_trans;
foreach (1 .. $param_N) {
  push (@p_trans, [@{$hash{"p_$_"}}]);
}

# Generate matrix of parameters p[m-1][n-1] and take their
# squareroot, thus allowing simplex parameters to be negative.
my @p;
for(my $i=0; $i<$ndim; $i++) {
  for(my $j=0; $j<$param_N; $j++) {
    $p[$i][$j]=sqrt($p_trans[$j][$i]);
  }
}

my @psum;
my @ptry;
foreach (1 .. $param_N) {
push(@ptry, sqrt(${$hash{"p_$_"}}[-1]));
}
my $ytry=$ftar[-1];

# Generate and sort arrays according to y[m-1]
my @i_sort;
my @y;
my @ftar_asc;

foreach (0 .. $param_N) {
  $y[$_]=$ftar[$_];
  $i_sort[$_]=$_;
}

@i_sort=(sort{$y[$a] <=> $y[$b]} @i_sort);

my @y_asc;
my @p_asc;
for (my $i=0;$i<=$#y;$i++) {
  $y_asc[$i]=$y[$i_sort[$i]];
  for (my $j=0;$j<$param_N;$j++){
    $p_asc[$i][$j]=$p[$i_sort[$i]][$j];
  }
}
@y=@y_asc;

# Define highest, next highest, and lowest points (ihi, inhi, ilo)
my $ihi=$#y;
my $ilo=0;
my $inhi=$ihi-1;
for (my $i=0;$i<=$#y;$i++) {
  if($y[$i]<=$y[$ilo]) {$ilo=$i;}
  if($y[$i]>$y[$ihi]) {
    $inhi=$ihi;
    $ihi=$i;
  }
  if ($y[$i]>$y[$inhi] && $i!=$ihi) {$inhi=$i;}
}

# Evalulate function at the trial point
if ($ytry<$y[$ihi]) {
  for (my $j=0;$j<$param_N;$j++) {
    $y[$ihi]=$ytry;
    $psum[$j]+=$ptry[$j]-$p[$ihi][$j];
    $p_asc[$ihi][$j]=$ptry[$j];
  }
}

# Get Transformation from previous state file
my %state;
open (STATE_CUR, "< state_$name.cur") || die "Could not open file state_$name.cur\n";
while(<STATE_CUR>) {
  if (/^(.*)=(.*)$/) {
  $state{"$1"}=$2;
  }
}
close(STATE_CUR);

# Prepare new state file
open (STATE, "> state_$name.new") || die "Could not open file state_$name.new\n";

# ---------------------------------------------------------------------------------
#                     DOWNHILL SIMPLEX ALGORITHM starting...                      |
# ---------------------------------------------------------------------------------

my $nfunc=$state{'nfunc'};

# If better than the best, try bigger step in the same direction (Expansion)
if ($state{'Transformation'} eq "Reflection" && $ytry <= $y[$ilo]) {
  print STATE "Transformation=Expansion\n";
  @psum=calc_psum(@p,$param_N,$ndim);
  @ptry=calc_ptry($param_N,$ihi,2.0,@p,@psum);
  push(@y,"0");
  my @empty=();
  push(@p_asc, \@empty);
  for (my $j=0;$j<$param_N;$j++){
    $p_asc[-1][$j]=$ptry[$j];
  }
  push(@flag,"pending");
  $nfunc++;
  print STATE "nfunc=$nfunc";
}

# if worse than the second worst, try smaller step (Contraction)
elsif ($state{'Transformation'} eq "Reflection" && $ytry >= $y[$inhi]) {
  print STATE "Transformation=Contraction\n";
  print STATE "ysave=$y[$ihi]\n";
  @psum=calc_psum(@p,$param_N,$ndim);
  @ptry=calc_ptry($param_N,$ihi,0.5,@p,@psum);
  push(@y,"0");
  my @empty=();
  push(@p_asc, \@empty);
  for (my $j=0;$j<$param_N;$j++){
    $p_asc[-1][$j]=$ptry[$j];
  }
push(@flag,"pending");
$nfunc++;
print STATE "nfunc=$nfunc";
}

# if worse than worst point from contraction, replace all but the best
elsif ($state{'Transformation'} eq "Contraction" && $ytry >= $state{'ysave'}) {
  print STATE "Transformation=Reduction\n";
  for (my $i=0;$i<$ndim;$i++) {
    if ($i!=$ilo) {
      for (my $j=0;$j<$param_N;$j++) {
        $p_asc[$i][$j]=$psum[$j]=0.5*($p[$i][$j]+$p[$ilo][$j]);
        $y[$i]="0";
        $flag[$i]="pending";
      }
    }
  }
  $mdim-=1;
  $nfunc+=$ndim;
  print STATE "nfunc=$nfunc";
}

else {
  # Compute reflected point
  print STATE "Transformation=Reflection\n";
  @psum=calc_psum(@p,$param_N,$ndim);
  @ptry=calc_ptry($param_N,$ihi,-1.0,@p,@psum);
  push(@y,"0");
  my @empty=();
  push(@p_asc, \@empty);
  for (my $j=0;$j<$param_N;$j++){
    $p_asc[-1][$j]=$ptry[$j];
  }
  push(@flag,"pending");
  $nfunc++;
  print STATE "nfunc=$nfunc";
}

# Check for convergence
my $rtol=2.0*abs($y[$ihi]-$y[$ilo])/(abs($y[$ihi])+abs($y[$ilo]));

if($rtol<$ftol) {
   ($y[$ilo],$y[0])=($y[0],$y[$ilo]);
      for (my $j=0;$j<$param_N;$j++){
      ($p[$ilo][$j],$p[0][$j])=($p[0][$j],$p[$ilo][$j]);
      print STATE "Done - Simplex converged after $nfunc steps.\n";
      die "--- Simplex convergerd after $nfunc steps ---";
   }
}

close(STATE);

for (my $i=0;$i<=$#y;$i++) {
  for (my $j=1;$j<=$param_N;$j++){
    ${$hash{"p_$j"}}[$i]=($p_asc[$i][$j-1])**2;
  }
}

# Update simplex table
saveto_simplex_table($outfile,$mdim,$param_N,@y,%hash,@flag) or die "$progname: error at saveto_simplex_table\n";
