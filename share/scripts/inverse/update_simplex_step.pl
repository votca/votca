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

# Generate matrix of parameters p[m-1][n-1] and take their
# squareroot, thus allowing simplex parameters to be negative.
my @p_trans;
my @p;
foreach (1 .. $param_N) {
  push (@p_trans, [@{$hash{"p_$_"}}]);
}
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

# Generate and arrays according to y[ilo]<...<y[inhi]<y[ihi]
my ($y_ref,$p_ref)=sort_ftar($param_N, @ftar, @p);
my @y_asc=@$y_ref;
my @p_asc=@$p_ref;

# Define highest, next highest, and lowest points (ihi, inhi, ilo)
my $ihi=$#y_asc;
my $ilo=0;
my $inhi=$ihi-1;

# Evalulate function at the trial point
if ($ytry<$y_asc[$ihi]) {
  for (my $j=0;$j<$param_N;$j++) {
    $y_asc[$ihi]=$ytry;
    $psum[$j]+=$ptry[$j]-$p_asc[$ihi][$j];
    $p_asc[$ihi][$j]=$ptry[$j];
  }
}

# Sort arrays according to y[ilo]<...<y[inhi]<y[ihi]
($y_ref,$p_ref)=sort_ftar($param_N, @y_asc, @p);
@y_asc=@$y_ref;
@p_asc=@$p_ref;

# Read previous state file
my %state_cur;
open (STATE_CUR, "< state_$name.cur") || die "Could not open file state_$name.cur\n";
while(<STATE_CUR>) {
  if (/^(.*)=(.*)$/) {
  $state_cur{"$1"}=$2;
  }
}
close(STATE_CUR);

# Prepare new state file
open (STATE_NEW, "> state_$name.new") || die "Could not create file state_$name.new\n";

# ---------------------------------------------------------------------------------
#                     DOWNHILL SIMPLEX ALGORITHM starting...                      |
# ---------------------------------------------------------------------------------

my $state_new;

switch ($state_cur{'Transformation'}) {

  case "None" {
    $state_new="Reflection";
  }
  case "Reflection" {
    if ($ytry <= $y_asc[$ilo]) {
      $state_new="Expansion";
    }
    elsif ($ytry >= $y_asc[$inhi]) {
      $state_new="Contraction";
    }
    else {
      $state_new="Reflection";
    }  
  }
  case "Expansion" {
    $state_new="Reflection";
  }
  case "Contraction" {
    if ($ytry >= $state_cur{'ysave'}) {
      $state_new="Reduction";
    }
    else {
        $state_new="Reflection";
    }
  }
  case "Reduction" {
    $state_new="Reflection";
  }
  else {
    die "error: undefined Transformation (state_new)\n";
  }

} # end of switch (state_cur)

my $nfunc=$state_cur{'nfunc'};

switch ($state_new) {

print STATE_NEW "Transformation=$state_new\n";

  case "Expansion" {
    @psum=calc_psum(@p_asc,$param_N,$ndim);
    @ptry=calc_ptry($param_N,$ihi,2.0,@p_asc,@psum);
    push(@y_asc,"0");
    my @empty=();
    push(@p_asc, \@empty);
    for (my $j=0;$j<$param_N;$j++){
      $p_asc[-1][$j]=$ptry[$j];
    }
    push(@flag,"pending");
    $nfunc++;
    print STATE_NEW "nfunc=$nfunc\n";
  }

  case "Contraction" {
    print STATE_NEW "ysave=$y_asc[$ihi]\n";
    @psum=calc_psum(@p_asc,$param_N,$ndim);
    @ptry=calc_ptry($param_N,$ihi,0.5,@p_asc,@psum);
    push(@y_asc,"0");
    my @empty=();
    push(@p_asc, \@empty);
    for (my $j=0;$j<$param_N;$j++){
      $p_asc[-1][$j]=$ptry[$j];
    }
    push(@flag,"pending");
    $nfunc++;
    print STATE_NEW "nfunc=$nfunc\n";
  }

  case "Reduction" {
    for (my $i=0;$i<$ndim;$i++) {
        if ($i!=$ilo) {
          for (my $j=0;$j<$param_N;$j++) {
            $p_asc[$i][$j]=$psum[$j]=0.5*($p_asc[$i][$j]+$p_asc[$ilo][$j]);
            $y_asc[$i]="0";
            $flag[$i]="pending";
          }
        }
      }
    $mdim-=1;
    $nfunc+=$ndim-1;
    print STATE_NEW "nfunc=$nfunc\n";
  }

  case "Reflection" {
    @psum=calc_psum(@p_asc,$param_N,$ndim);
    @ptry=calc_ptry($param_N,$ihi,-1.0,@p_asc,@psum);
    push(@y_asc,"0");
    my @empty=();
    push(@p_asc, \@empty);
    for (my $j=0;$j<$param_N;$j++){
      $p_asc[-1][$j]=$ptry[$j];
    }
    push(@flag,"pending");
    $nfunc++;
    print STATE_NEW "nfunc=$nfunc\n";
  }

} # end of switch (state_new)

# Check for convergence
my $rtol=2.0*abs($y_asc[$ihi]-$y_asc[$ilo])/(abs($y_asc[$ihi])+abs($y_asc[$ilo]));

if($rtol<$ftol) {
   ($y_asc[$ilo],$y_asc[0])=($y_asc[0],$y_asc[$ilo]);
      for (my $j=0;$j<$param_N;$j++){
      ($p_asc[$ilo][$j],$p_asc[0][$j])=($p_asc[0][$j],$p_asc[$ilo][$j]);
      die "--- Simplex convergerd after $nfunc steps ---";
   }
}

close(STATE_NEW);

for (my $i=0;$i<=$#y_asc;$i++) {
  for (my $j=1;$j<=$param_N;$j++){
    ${$hash{"p_$j"}}[$i]=($p_asc[$i][$j-1])**2;
  }
}

# Update simplex table
saveto_simplex_table($outfile,$mdim,$param_N,@y_asc,%hash,@flag) or die "$progname: error at saveto_simplex_table\n";
