#! /usr/bin/perl -w

use strict;

( my $progname = $0 ) =~ s#^.*/##;

if ("$ARGV[0]" eq "--help"){
  print <<EOF;
Usage: $progname infile1 infile2 outfile
This script adds up two potentils 
In addtion it does some magic tricks:
 - order of infiles MATTER !!!!
 - if infile2 contain undef value, it uses the value from infile1
 - if value for infile1 and infile2 are invalid, result is also invalid
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
my @flag;

# TODO: think about addition rules
# now I did it like that to always maintain interval of interest in all potentials
# shount that just be a < instead of <= ??
for (my $i=0;$i<=$#r_cur;$i++){
  if($flag_cur[$i] eq "u" || $flag_delta[$i] eq "u") {
    $pot[$i] = $pot_cur[$i];  # is already nan or we don't change
    $flag[$i] = "u";
  }
  else {
    $pot[$i]=$pot_cur[$i]+$pot_delta[$i];
    $flag[$i] = $flag_cur[$i];
  }
  #if ($flag_cur[$i] eq "i"){
  #  if ($flag_delta[$i] eq "i"){
  #    $pot[$i]=$pot_cur[$i]+$pot_delta[$i];
  #  } else {
  #    $pot[$i]=$pot_cur[$i];
  #  }
  #  $flag[$i]="i";
  #} else {
  #  $pot[$i]="nan";
  #  $flag[$i]="u";
  #}
}
saveto_table($outfile,\@r_cur,\@pot,\@flag) || die "$progname: error at save table\n";

