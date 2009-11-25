#! /usr/bin/perl -w
use strict;

( my $progname = $0 ) =~ s#^.*/##;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
Usage: $progname infile outfile
This script set the beginning of the dpot to the first valid value and shift the whole potential
so that dpot(r_max)=0.

NEEDS:
USES: readin_table saveto_table 
EOF
  exit 0;
}

die "2 parameters are nessary, <infile> <outfile>\n" if ($#ARGV<1);

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

print("$infile, $outfile\n");

# read in the current dpot
my @r;
my @dpot;
my @flag;
(readin_table($infile,@r,@dpot,@flag)) || die "$progname: error at readin_table\n";

# find first u/o
my $i_first;
for($i_first=0; ($i_first<=$#r) && ($flag[$i_first] =~ /[uo]/); $i_first++) {}

# shift beginning
for(my $i=0; $i<$i_first; $i++) {
    $dpot[$i] = $dpot[$i_first];
    $flag[$i]="o"
}

# bring end to zero
for(my $i=0; $i<=$#r; $i++) {
    $dpot[$i] -= $dpot[$#r];
}

# save to file
saveto_table($outfile,@r,@dpot,@flag) || die "$progname: error at save table\n";
