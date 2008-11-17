#! /usr/bin/perl -w

use strict;

die "2 parameters are nessary\n" if ($#ARGV<1);

my $file="$ARGV[0]";
open(FILE1,$file) or die "$file not found\n";

my $file2="$ARGV[1]";
open(FILE2,"> $file2") or die "$file2 not found\n";

my @r;
my @pot;
while (<FILE1>){
   next if /^#/;
   my @parts=split(' ');
   push(@r,$parts[0]);
   push(@pot,$parts[1]);
}

my $gro_delta=0.002;

my $tmp_r=0;
my $r_cut=0.9;
my $r_max=2;

while ($tmp_r<$r[0]){
   print FILE2 "$tmp_r 0.0 0.0 0.0 0.0 0.0 0.0\n";
   $tmp_r+=$gro_delta;
}

#first line no derivative
print FILE2 "$r[0] 0.0 0.0 0.0 0.0 $pot[0] 0.0\n";

my $deri_pot;
for (my $i=1;$i<$#r;$i++){
   $deri_pot=-($pot[$i+1]-$pot[$i-1])/($r[$i+1]-$r[$i-1]);
   print FILE2 "$r[$i] 0.0 0.0 0.0 0.0 $pot[$i] $deri_pot\n";
}

#last point
print FILE2 "$r[$#pot] 0.0 0.0 0.0 0.0 $pot[$#pot] 0.0\n";

#add zeros
$tmp_r=$r[$#pot]+$gro_delta;
while($tmp_r<$r_max+$gro_delta){
   print FILE2 "$tmp_r 0.0 0.0 0.0 0.0 0.0 0.0\n";
   $tmp_r+=$gro_delta;
}
