#! /usr/bin/perl -w

use strict;

( my $progname = $0 ) =~ s#^.*/##;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
Usage: $progname infile outfile
This script convert csg potential files to xvg format
Potential are copy in the C12 column 
In addtion it does some magic tricks:
 -bigger value will be set to pot_max (see xml)
 -shift the potential, so that it is zero at the cutoff
 -set all values to zero after the cutoff

NEEDS: cg.inverse.gromacs.pot_max cg.inverse.gromacs.table_end cg.inverse.gromacs.table_bins
USES: \$SOURCE_WRAPPER csg_get_property saveto_table readin_table
EOF
  exit 0;
}

(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";

die "2 parameters are nessary\n" if ($#ARGV<1);

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my $gromacs_max=csg_get_property("cg.inverse.gromacs.pot_max");
my $table_end=csg_get_property("cg.inverse.gromacs.table_end");
my $table_bins=csg_get_property("cg.inverse.gromacs.table_bins");

my @r;
my @pot;
my @flag;
(readin_table($infile,\@r,\@pot,\@flag)) || die "$progname: error at readin_table\n";

#gromacs does not like VERY big numbers
for (my $i=0;$i<=$#r;$i++) {
  $pot[$i]=$gromacs_max if $pot[$i]>$gromacs_max;
}

#cutoff is last point
my $i_cut=$#r;

#shift potential so that it is zero at cutoff
for (my $i=0;$i<=$i_cut;$i++){
   $pot[$i]-=$pot[$i_cut];
}

# set end of the potential to zero
for (my $i=$i_cut;$i<=$table_end/$table_bins;$i++) {
  $pot[$i]=0;
  $r[$i]=$r[$i-1]+$table_bins;
}

my @force;

#calc force
$force[0]=0;
for (my $i=1;$i<$#r;$i++){
   $force[$i]=-($pot[$i+1]-$pot[$i-1])/($r[$i+1]-$r[$i-1]);
}
$force[$#r]=0.0;

open(OUTFILE,"> $outfile") or die "saveto_table: could not open $outfile\n";
for(my $i=0;$i<=$#r;$i++){
  print OUTFILE "$r[$i] 0.0 0.0 0.0 0.0 $pot[$i] $force[$i]\n";
}
close(OUTFILE) or die "Error at closing $outfile\n";

