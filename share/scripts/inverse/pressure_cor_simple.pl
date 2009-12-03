#! /usr/bin/perl -w

use strict;
( my $progname = $0 ) =~ s#^.*/##;
if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
Usage: $progname p_cur outfile
This script calls the pressure corrections dU=A*(1-r/r_c), where
A=-0.1k_B T * max(1,|p_cur-p_target|*scale) * sign(p_cur-p_target)

NEEDS: cg.inverse.kBT max step inverse.p_target inverse.post_update_options.pressure.simple.scale
USES: csg_get_property csg_get_interaction_property saveto_table
EOF
  exit 0;
}

die "2 parameters are nessary\n" if ($#ARGV<1);

use CsgFunctions;

my $kBT=csg_get_property("cg.inverse.kBT");
my $max=csg_get_interaction_property("max");
my $delta_r=csg_get_interaction_property("step");
my $scale_factor=csg_get_interaction_property("inverse.post_update_options.pressure.simple.scale");
my $p_target=csg_get_interaction_property("inverse.p_target");
my $p_now=$ARGV[0];

#Determine the sign
my $pref;
if ($p_now>$p_target){
   $pref=-0.1*$kBT;
} else {
   $pref=0.1*$kBT;
}

#Determine pressure factor
my $p_factor=($p_now-$p_target)*$scale_factor;
$p_factor=-$p_factor if $p_factor<0;

#Only use pressure factor if not too big
#max is 0.1kbT
$pref*=$p_factor if $p_factor<1;

my @r;
my @pot;
my @flag;
my $outfile="$ARGV[1]";
for(my $i=0;$i<=$max/$delta_r;$i++){
  $r[$i]=$i*$delta_r;
  $pot[$i]=$pref*(1-$r[$i]/$max);
  $flag[$i]="i";
}
saveto_table($outfile,@r,@pot,@flag) || die "$progname: error at save table\n";

