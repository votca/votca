#! /usr/bin/perl -w

use strict;


# added this hope its righ (victor)
( my $progname = $0 ) =~ s#^.*/##;
(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";
####

die "3 parameters are nessary\n" if ($#ARGV<2);

my $pref=get_sim_property("kBT");
my $r_cut=csg_get("cut");
my $r_max=csg_get("max");
my $delta_r=csg_get("step");

#my $pref=300*0.00831451;
#my $r_cut=0.9;
#my $delta_r=0.01;

my $small_x=0.001;

my $file="$ARGV[0]";
open(FILE1,$file) or die "$file not found\n";
my @r_aim;
my @rdf_aim;
while (<FILE1>){
   next if /^[#@]/;
   my @parts=split(' ');
   push(@r_aim,$parts[0]);
   push(@rdf_aim,$parts[1]);
}
close(FILE1) or die "Error at closing $file\n";

$file="$ARGV[1]";
open(FILE1,$file) or die "$file not found\n";
my @r_cur;
my @rdf_cur;
while (<FILE1>){
   next if /^[#@]/;
   my @parts=split(' ');
   push(@r_cur,$parts[0]);
   push(@rdf_cur,$parts[1]);
}
close(FILE1) or die "Error at closing $file\n";

die "Different grids \n" if (($r_aim[1]-$r_aim[0])!=($r_cur[1]-$r_cur[0]));

die "Different start point \n" if (($r_aim[0]-$r_cur[0]) > 0.0);

#find min_i_aim
my $i_min_aim;
for (my $i=$#r_aim;$i>=0;$i--){
   $i_min_aim=$i+1;
   last if $rdf_aim[$i]<$small_x;
}

#find min_i_cut
my $i_min_cur;
for (my $i=$#r_cur;$i>=0;$i--){
   $i_min_cur=$i+1;
   last if $rdf_cur[$i]<$small_x;
}

my $i_min=($i_min_aim>$i_min_cur)?$i_min_aim:$i_min_cur;

#print "XXX $i_min\n";
$file="$ARGV[2]";
open(FILE,"> $file") or die "$file not found\n";

for (my $nr=$i_min;$nr<=$#r_cur;$nr++){
   my $tmp_pot;
   if($rdf_aim[$nr] != 0) {
      $tmp_pot=log($rdf_cur[$nr]/$rdf_aim[$nr])*$pref;
   }
   else {
    $tmp_pot=0;
   }

   print FILE "$r_cur[$nr] $tmp_pot\n";
   #print "YYY $nr $tmp_r $tmp_pot\n";
}

close(FILE);
