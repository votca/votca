#! /usr/bin/perl -w

use strict;

( my $progname = $0 ) =~ s#^.*/##;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
Usage: $progname infile outfile
This script smoothes a table

USES: readin_table saveto_table
NEEDS: 
EOF
  exit 0;
}

die "2 parameters are nessary\n" if ($#ARGV<1);

use CsgFunctions;
(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";

my $infile="$ARGV[0]";
my @r_cur;
my @pot_cur;
my @flag_cur;
(readin_table($infile,@r_cur,@pot_cur,@flag_cur)) || die "$progname: error at readin_table\n";

my $outfile="$ARGV[1]";
my @pot;

# TODO: think about addition rules
# now I did it like that to always maintain interval of interest in all potentials
for (my $i=1;$i<$#r_cur;$i++){
  $pot[$i]=$pot_cur[$i]; 
  if($flag_cur[$i] eq "i") {
    $pot[$i] = 0.25*$pot_cur[$i-1] + 0.5*$pot_cur[$i] + 0.25*$pot_cur[$i+1]; 
  }
}

$pot[0]=$pot_cur[0];
$pot[$#pot_cur]=$pot_cur[$#pot_cur];

if($flag_cur[0] eq "i") {
  $pot[0] = (2.*$pot_cur[0] + $pot_cur[1])/3;
}
if($flag_cur[$#pot_cur] eq "i") {
  $pot[$#pot_cur] = (2.*$pot_cur[$#pot_cur] + $pot_cur[$#pot_cur-1])/3;
}

saveto_table($outfile,@r_cur,@pot,@flag_cur) || die "$progname: error at save table\n";

