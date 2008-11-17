#! /usr/bin/perl -w

use strict;

die "3 parameters are nessary\n" if ($#ARGV<2);

my $file="$ARGV[0]";
open(FILE1,$file) or die "$file not found\n";
my @r_cur;
my @pot_cur;
while (<FILE1>){
   next if /^#/;
   my @parts=split(' ');
   push(@r_cur,$parts[0]);
   push(@pot_cur,$parts[1]);
}
close(FILE1) or die "Error at closing $file\n";

$file="$ARGV[1]";
open(FILE1,$file) or die "$file not found\n";
my @r_delta;
my @pot_delta;
while (<FILE1>){
   next if /^#/;
   my @parts=split(' ');
   push(@r_delta,$parts[0]);
   push(@pot_delta,$parts[1]);
}
close(FILE1) or die "Error at closing $file\n";

die "Different grids $r_delta[1] $r_delta[0] $r_cur[1] $r_cur[0]\n" if (($r_delta[1]-$r_delta[0]-$r_cur[1]+$r_cur[0])>0.0001);

my $delta=$r_delta[1]-$r_delta[0];
my $shift=int(($r_delta[0]-$r_cur[0])/$delta)+1;
#print "XXX $delta $shift $r_cur[$shift] $r_delta[0]\n";

my $pot_shift=$pot_delta[$#r_cur-$shift]*exp(-1);

print "XXX $pot_shift\n";
$file="$ARGV[2]";
open(FILE,"> $file") or die "$file not found\n";
for (my $i=0;$i<=$#r_cur;$i++){
   my $tmp;
   if ($i>=$shift){
      $tmp=$pot_cur[$i]+$pot_delta[$i-$shift];
      #print "xxx $i $r_cur[$i] $pot_cur[$i] $pot_delta[$i-$shift] $tmp\n";
   } else {
      $tmp=$pot_cur[$i];
   }
   $tmp*=exp(-($i-$#r_cur+5)/5.) if ($i>$#r_cur-5);
   #printf "xxx $i %f\n",exp(-($i-$#r_cur+5)) if ($i>=$#r_cur-5);
   $tmp-=$pot_shift;
   print FILE "$r_cur[$i] $tmp\n";
}
close(FILE);
