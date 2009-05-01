#! /usr/bin/perl -w

use strict;

( my $progname = $0 ) =~ s#^.*/##;

if ("$ARGV[0]" eq "--help"){
  print <<EOF;
Usage: $progname infile1 infile2 outfile
This script adds up two potentils 
In addtion it does some magic tricks:
 -do not crash when calc log(0)
 -extrapolate the beginnig of pot
 -the maximum to interpolate is pot_max (see xml)
 -bigger value will be set to that max
 -shift the potential, so that it is zero at the cutoff
 -set all values to zero after the cutoff
EOF
  exit 0;
}

(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";

die "3 parameters are nessary\n" if ($#ARGV<2);

my $infile="$ARGV[0]";
my @r_cur;
my @pot_cur;
my @flag_cur;
(readin_table($infile,\@r_cur,\@pot_cur,\@flag_cur)) || die "$progname: error at readin_table\n";

my $infile2="$ARGV[1]";
#delta is just a name
my @r_delta;
my @pot_delta;
my @flag_delta;
(readin_table($infile2,\@r_delta,\@pot_delta,\@flag_delta)) || die "$progname: error at readin_table\n";

#should never happen, but ....
die "Different grids\n" if (($r_delta[1]-$r_delta[0]-$r_cur[1]+$r_cur[0])>0.0001);
die "Different start point \n" if (($r_delta[0]-$r_cur[0]) > 0.0);

my $outfile="$ARGV[2]";
my @pot;
for (my $i=0;$i<=$#r_cur;$i++){
  if (($flag_delta[$i] eq "i" ) && ($flag_cur[$i] eq "i")){
    $pot[$i]=$pot_cur[$i]+$pot_delta[$i];
  } else {
    $pot[$i]=0;
    $flag_cur[$i]="u";
  }
}
saveto_table($outfile,\@r_cur,\@pot,\@flag_cur) || die "$progname: error at save table\n";

