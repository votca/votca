#! /usr/bin/perl -w
#
# shifts the beginning of a non-bonded delta potential which is not in the update
# region to same as first value which is in region. It also shits the whole dpot 
# that it is zero at r_cut
#
# TODO: care about flags, not tested
#
use strict;

# source the function file
( my $progname = $0 ) =~ s#^.*/##;
(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";
######

die "2 parameters are nessary, <infile> <outfile>\n" if ($#ARGV<2);

my $file="$ARGV[0]";
open(FILE1,$file) or die "$file not found\n";

my $r_cut=csg_get("max");

# read in the current dpot
my @r;
my @dpot;
#my @flag;
while (<FILE1>){
   next if /^#/;
   my @parts=split(' ');
   push(@r,$parts[0]);
   push(@dpot,$parts[1]);
#   push(@flag,$parts[2]);
}
close(FILE1) or die "Error at closing $file\n";

# find index of cutoff
my $i_cut;
for($i_cut=0; ($i_cut<$#r) && ($r[$i_cut]<$r_cut); $i_cut++) {}
# find last u/o
my $i_first;
for(my $i_first=0; ($i<$#r) && ($dpot[$i_first] == 0); $i_first++) {}

# shift beginning
for(my $i=0; $i<$i_first; $i++) {
    $dpot[$i] = $dpot[$i_first];
}

# bring end to zero
for(my $i=0; $i<$i_first; $i++) {
    $dpot[$i] -= $dpot[$i_cut];
}

# save to file
$file="$ARGV[2]";
open(FILE,"> $file") or die "$file not found\n";
for (my $i=0;$i<=$#r_cur;$i++) {
   print FILE "$r[$i] $dpot[$i] i\n";
}

close(FILE);
