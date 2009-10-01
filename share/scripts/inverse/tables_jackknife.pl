#! /usr/bin/perl -w

use strict;

( my $progname = $0 ) =~ s#^.*/##;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
Usage: $progname out full block1 block2 ...
This scripts calculates the jacknife error from existing tables 
full: table calculated with full dataset
blocks: tables calculated with 1 block missing
outfile: file to write results

USES: \$SOURCE_WRAPPER readin_table saveto_table_err
NEEDS: 
EOF
  exit 0;
}

(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";

die "3 parameters are nessary\n" if ($#ARGV<2);

my $file_full="$ARGV[1]";
my @r_full;
my @val_full;
my @flag_full;

my $outfile="$ARGV[0]";
my @err;


(readin_table($file_full,\@r_full,\@val_full,\@flag_full)) || die "$progname: error at readin_table\n";

for (my $i=0;$i<=$#r_full;$i++) {
  $err[$i]=0;
}

shift @ARGV;
shift @ARGV;

my $nblocks = 0;
while (@ARGV > 0) {
  my $file_cur="$ARGV[0]";
  my @r_cur;
  my @val_cur;
  my @flag_cur;

  (readin_table($file_cur,\@r_cur,\@val_cur,\@flag_cur)) || die "$progname: error at readin_table\n";
  #should never happen, but ....
  #die "Different grids\n" if (($r_delta[1]-$r_delta[0]-$r_cur[1]+$r_cur[0])>0.0001);
  #die "Different start point \n" if (($r_delta[0]-$r_cur[0]) > 0.0);
 
  for (my $i=0;$i<=$#r_cur;$i++) {
      $err[$i] += ($val_cur[$i] - $val_full[$i])**2;  # is already nan or we don't change
  }
  shift @ARGV;
  $nblocks = $nblocks + 1;
}

for (my $i=0;$i<=$#r_full;$i++) {
  $err[$i]=sqrt(($nblocks-1)/$nblocks*$err[$i]);
}

saveto_table_err($outfile,\@r_full,\@val_full,\@flag_full,\@err) || die "$progname: error at save table\n";

