#! /usr/bin/perl -w

use strict;

( my $progname = $0 ) =~ s#^.*/##;
(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";

die "2 parameters are nessary\n" if ($#ARGV<1);

my $file="$ARGV[0]";
open(FILE1,$file) or die "$file not found\n";

my $file2="$ARGV[1]";
open(FILE2,"> $file2") or die "$file2 not found\n";

my $pref=get_sim_property("kBT");
my $r_cut=csg_get("cut");
# my $r_max=csg_get("max"); # these are not needed (vr)
# my $delta_r=csg_get("step"); # these are not needed (vr)

my @r;
my @pot;
my @flag;

while (<FILE1>){
   next if /^[#@]/;
   my @parts=split(' ');
   push(@r,$parts[0]);
   my $tmp=$parts[1]>0.0?-$pref*log($parts[1]):"0";
   push(@pot,$tmp);
   my $tmp_fl=$parts[1]>0.0?"i":"u";
   push(@flag,$tmp_fl);
}

#find last nan
my $nr;
for ($nr=$#pot-1;$nr>=0;$nr--){
   last if ($flag[$nr] =~ /u/) ;
}
$nr++;
my $pot_max=$pot[$nr];

#extra polation the begin
my $deri=$pot[$nr]-$pot[$nr+1];
while ($nr>=0) { #$pot_max<100000){
   $deri+=$deri;
   $pot_max=$pot[$nr]+$deri;
   $nr--;
   $pot[$nr]=$pot_max<100000?$pot_max:100000;
   $flag[$nr]=$pot_max<100000?"o":"u";
}
my $i_start=$nr;

#find r_cut
for ($nr=0;$nr<=$#r;$nr++){
   last if $r[$nr]>=$r_cut;
}
#print "XXX $pot[$nr] $nr\n";
my $i_cut=$nr;

#substract potential of r_cut
for (my $i=0;$i<=$i_cut;$i++){
   $pot[$i]-=$pot[$i_cut] unless  ($flag[$i] =~ /[uo]/);
}

#smooth end
for (my $i=$i_cut-5;$i<=$i_cut;$i++){
   #print "YYY $i $r[$i] $pot[$i]\n";
   $pot[$i]=$pot[$i]*exp(-($r[$i]/$r[$i_cut-5]-1)**2);
   #print "ZZZ $i $r[$i]  $pot[$i]\n";
}

# set ends to zero
for (my $i=$i_cut;$i<$#flag;$i++) {
  $pot[$i]=0;
  $flag[$i]="o";
}

for(my $count=0;$count<$#pot;$count++){
   print FILE2 "$r[$count] $pot[$count] $flag[$count]\n";
}

close(FILE1) or die "Error at closing $file\n";
close(FILE2) or die "Error at closing $file2\n";
