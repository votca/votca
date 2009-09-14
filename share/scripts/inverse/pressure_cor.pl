#! /usr/bin/perl -w

use strict;
( my $progname = $0 ) =~ s#^.*/##;
if ("$ARGV[0]" eq "--help"){
  print <<EOF;
Usage: $progname p_target p_cur outfile
This script calls the pressure corrections dU=A*(1-r/r_c) 
EOF
  exit 0;
}

(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";

die "3 parameters are nessary\n" if ($#ARGV<2);

my $kBT=csg_get_property("cg.inverse.kBT");
my $max=csg_get("max");
my $delta_r=csg_get("step");

my $p_target=$ARGV[0];
my $p_now=$ARGV[1];

#Determine the sign
my $pref;
if ($p_now>$p_target){
   $pref=-0.1*$kBT;
} else {
   $pref=0.1*$kBT;
}

#Determine pressure factor
my $p_factor=($p_now-$p_target)/3000;
$p_factor=-$p_factor if $p_factor<0;

#Only use pressure factor if not too big
#max is 0.1kbT
$pref*=$p_factor if $p_factor<1;

my @r;
my @pot;
my @flag;
my $outfile="$ARGV[2]";
for(my $i=0;$i<=$max/$delta_r;$i++){
  $r[$i]=$i*$delta_r;
  $pot[$i]=$pref*(1-$r[$i]/$max);
  $flag[$i]="i";
}
saveto_table($outfile,\@r,\@pot,\@flag) || die "$progname: error at save table\n";

