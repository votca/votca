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
# ---------------------------
# Downhill Simplex Algorithm 
# ---------------------------

use strict;

( my $progname = $0 ) =~ s#^.*/##;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
   print <<EOF;
$progname, version %version%
This script adds up two potentils 
In addtion it does some magic tricks:
- order of infiles MATTER !!!!
- if infile2 contain undef value, it uses the value from infile1
- if value for infile1 and infile2 are invalid, result is also invalid

Usage: $progname infile1 infile2 outfile

NEEDS:

USES: readin_table saveto_table 
EOF
  exit 0;
}

die "2 parameters are nessary\n" if ($#ARGV<1);

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

use CsgFunctions;
use SimplexFunctions;
use Switch;

my @ftar;
my @sig;
my @eps;
my @flag_simplex;
(readin_simplex_table($infile,@ftar,@sig,@eps,@flag_simplex)) || die "error at readin_simplex_table\n";

my $mpts=3;
my $ndim=$mpts-1;
my $nfunc=0;
my $ftol=csg_get_property("cg.inverse.ftol");

my @psum;
my @ptry;
my $ytry=$ftar[-1];

# Generate p[mpts][ndim] matrix (parameters)
my @p;
my @p_trans=([@sig],[@eps]);
for(my $i=0; $i<$mpts; $i++) {
   for(my $j=0; $j<$ndim; $j++) {
      $p[$i][$j]=$p_trans[$j][$i];
   }
}

# Generate and sort y[mpts] array (ftar values)
my $ilo=0;
my @i_sort;
my @ftar_asc;
my @sig_asc;
my @eps_asc;
foreach (0 .. $#ftar) {$i_sort[$_]=$_};
   @i_sort=(sort{$ftar[$a] <=> $ftar[$b]} @i_sort);
      for (my $i=0;$i<$mpts;$i++) {
      $ftar_asc[$i]=$ftar[$i_sort[$i]];
      $sig_asc[$i]=$sig[$i_sort[$i]];
      $eps_asc[$i]=$eps[$i_sort[$i]];
      }

# Define highest, next highest, and lowest points (ihi, inhi, ilo)
my @y=@ftar_asc;
my $ihi=$#ftar_asc;
my $inhi=0;
for (my $i=0;$i<$mpts;$i++) {
   if($y[$i]<=$y[$ilo]) {$ilo=$i;}
   if($y[$i]>$y[$ihi]) {
        $inhi=$ihi;
        $ihi=$i;
        }
   if ($y[$i]>$y[$inhi] && $i!=$ihi) {$inhi=$i;}
}

my %state;
open (STATE_CUR, "<state.cur") || die "Could not open file $_[0]\n";
while(<STATE_CUR>) {
   # create a hash
   if (/^(.*)=(.*)$/) {
   $state{"$1"}=$2;
   }
}
close(STATE_CUR);

open (STATE, ">state.new") || die "Could not open file $_[0]\n";

switch ($state{'Transformation'}) {

my $ysave_R;
my $ysave_C;
case 'Reflection' {
   # If better than the best, try bigger step in the same direction (Expansion)
   if ($ytry <= $y[$ilo]) {
      $ysave_R=$state{$ytry};
      print STATE "Transformation=Expansion\n";
      @psum=calc_psum(@p,$mpts,$ndim);
      @ptry=calc_ptry($ndim,$ihi,2.0,@p,@psum);
      push(@ftar_asc,"0");
      push(@sig_asc,"$ptry[0]");
      push(@eps_asc,"$ptry[1]");
      push(@flag_simplex,"pending");
      $nfunc++;
   }
   # if worse than the second worst as well as the worst, try smaller step (Contraction)
   elsif ($ytry >= $y[$inhi] && $ytry >= $y[$ihi]) {
      $ysave_C=$state{$y[$ihi]};
      print STATE "Transformation=Contraction\n";
      @psum=calc_psum(@p,$mpts,$ndim);
      @ptry=calc_ptry($ndim,$ihi,0.5,@p,@psum);
      push(@ftar_asc,"0");
      push(@sig_asc,"$ptry[0]");
      push(@eps_asc,"$ptry[1]");
      push(@flag_simplex,"pending");
      $nfunc++;
   }
   else {
   # Otherwise, replace worst point by reflected point
         for (my $j=0;$j<$ndim;$j++) {
         $y[$ihi]=$ytry;
         $p[$j]+=$ptry[$j]-$p[$ihi][$j];
         $p[$ihi][$j]=$ptry[$j];
      }
   }
}

case 'Expansion' {
   # If better the best, replace worst point by expanded point
   if ($ytry <= $y[$ilo]) {
      for (my $j=0;$j<$ndim;$j++) {
         $y[$ihi]=$ytry;
         $p[$j]+=$ptry[$j]-$p[$ihi][$j];
         $p[$ihi][$j]=$ptry[$j];
      }
   }
   else {
   # Otherwise, replace worst point by reflected point
      for (my $j=0;$j<$ndim;$j++) {
         $y[$ihi]=$ysave_R;
         $p[$j]+=$ptry[$j]-$p[$ihi][$j];
         $p[$ihi][$j]=$ptry[$j];
      }
   }
}

case 'Contraction' {
   # if worse than worst point from contraction, replace all but the best 
   if ($ytry>=$ysave_C) {
      print STATE "Transformation=Reduction\n";
      for (my $i=0;$i<$mpts;$i++) {
         if ($i!=$ilo) {
            for (my $j=0;$j<=$ndim;$j++) {
            $p[$i][$j]=0.5*($p[$i][$j]+$p[$ilo][$j]);
            $ftar_asc[$i]="0";
            $sig_asc[$i]=$p[$i][0];
            $eps_asc[$i]=$p[$i][1];
            $flag_simplex[$i]="pending";
            }
         }
      }
   $nfunc+=$ndim;
   }
   # Otherwise, replace worst point by contracted point
   else { 
      for (my $j=0;$j<$ndim;$j++) {
         $y[$ihi]=$ytry;
         $p[$j]+=$ptry[$j]-$p[$ihi][$j];
         $p[$ihi][$j]=$ptry[$j];
      }
   }
}

else { 
   # Compute reflected point
   print STATE "Transformation=Reflection\n";
   @psum=calc_psum(@p,$mpts,$ndim);
   @ptry=calc_ptry($ndim,$ihi,-1.0,@p,@psum);
   push(@ftar_asc,"0");
   push(@sig_asc,"$ptry[0]");
   push(@eps_asc,"$ptry[1]");
   push(@flag_simplex,"pending");
   $nfunc++;
}

} # End of switch loop

# Check for convergence
my $rtol=2.0*abs($y[$ihi]-$y[$ilo])/(abs($y[$ihi])+abs($y[$ilo]));

if($rtol<$ftol) {
   ($y[$ilo],$y[0])=($y[0],$y[$ilo]);
      for (my $i=0;$i<$ndim;$i++){
      ($p[$ilo][$i],$p[0][$i])=($p[0][$i],$p[$ilo][$i]);
      print STATE "Done - Simplex converged after $nfunc steps.\n";
      die "--- Simplex convergerd after $nfunc steps ---";
   }
}

close(STATE);

# Update simplex table
saveto_simplex_table($outfile,@ftar_asc,@sig_asc,@eps_asc,@flag_simplex) || die "$progname: error at save table\n";
