#! /usr/bin/perl -w

use strict;

# added this hope its righ (victor)
( my $progname = $0 ) =~ s#^.*/##;
(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";
####

die "3 parameters are nessary\n" if ($#ARGV<2);

my $kBT=get_sim_property("kBT");
my $max=csg_get("max");
my $delta_r=csg_get("step");

my $p_target=$ARGV[0];
my $p_now=$ARGV[1];

# my $file="rdf.xvg";
# open(FILE1,$file) or die "$file not found\n";
# my $integral=0;
# while (<FILE1>){
#    next if /^#/;
#    my @parts=split(' ');
#    my $rdf=$parts[1];
#    my $r=$parts[0];
#    if (($rdf>=0) && ($r<=$r_cut)){
#       $integral+=$rdf*$r*$r*$r*$delta_r;
#    }
# }
# close(FILE1) or die "Error at closing $file\n";


my $pref;
if ($p_now>$p_target){
   $pref=-0.1*$kBT;
} else {
   $pref=0.1*$kBT;
}

my $p_factor=($p_now-$p_target)/3000;
$p_factor=-$p_factor if $p_factor<0;
$pref*=$p_factor if $p_factor<1;

my $file="$ARGV[2]";
open(FILE,"> $file") or die "$file not found\n";
for(my $count=0;$count<=$max/$delta_r;$count++){
   my $tmp_r=$count*$delta_r;
   my $pot=$pref*(1-$tmp_r/$max);
   print FILE "$tmp_r $pot\n";
}
close(FILE);
