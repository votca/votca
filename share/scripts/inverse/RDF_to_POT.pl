#! /usr/bin/perl -w

use strict;

( my $progname = $0 ) =~ s#^.*/##;

if ("$ARGV[0]" eq "--help"){
  print <<EOF;
Usage: $progname infile outfile
This script convert in rdf to pot of mean force (F(r)=-k_B T*ln(g(r))
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

die "2 parameters are nessary\n" if ($#ARGV<1);

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

# TODO: this gromacs option should not be here 
#       since it's a general initial guess files
#       move this option out of gromacs section!!!!!
#
my $gromacs_max=csg_get_property("cg.inverse.gromacs.pot_max");
my $pref=csg_get_property("cg.inverse.kBT");
my $r_cut=csg_get("max");

my @r;
my @rdf;
my @flag;
(readin_table($infile,\@r,\@rdf,\@flag)) || die "$progname: error at readin_table\n";

my @pot;
for (my $i=0;$i<=$#r;$i++){
#  if ($flag[$i] eq "i"){
    #rdf = 0 will give undefined pot 
    if ($rdf[$i]>1e-10) {
      $pot[$i]=-$pref*log($rdf[$i]);
    }
    else {
      $pot[$i]="nan";
      $flag[$i]="u";
    }
#  }
}

#find first defined value (begining for r=0)
#but it is more stable to search first undefined value begin
#beginning form large r
my $first_undef_bin;
for (my $i=$#pot;$i>=0;$i--){
   if ($flag[$i] eq "u") {
     $first_undef_bin=$i;
     last;
   }
}

#gromacs does not like VERY big numbers
#in the very rare case that we are already in this region
#we try to find a new beginnig
while ($pot[$first_undef_bin+1]>$gromacs_max){
  $pot[$first_undef_bin+1]="nan";
  $flag[$first_undef_bin+1]="u";
  $first_undef_bin++;
}

#find i which is the cutoff
my $i_cut;
for (my $nr=0;$nr<=$#r;$nr++){
   if ($r[$nr]>=$r_cut) {
     $i_cut=$nr;
     last;
   }
}

#shift potential so that it is zero at cutoff
#first do the shift, then extrapolation
for (my $i=0;$i<=$i_cut;$i++){
   $pot[$i]-=$pot[$i_cut] unless  ($flag[$i] =~ /[u]/);
}

#quadratic extrapolation at the begining
#and set all undef values to max
my $slope=$pot[$first_undef_bin+1]-$pot[$first_undef_bin+2];
for (my $i=$first_undef_bin;$i>=0;$i--){
   $slope+=$slope;
   $pot[$i]=($pot[$i+1]+$slope)>$gromacs_max?$gromacs_max:($pot[$i+1]+$slope);
   $flag[$i]="o";
}

# set end of the potential to zero
for (my $i=$i_cut;$i<$#flag;$i++) {
  $pot[$i]=0;
  $flag[$i]="o";
}

saveto_table($outfile,\@r,\@pot,\@flag) || die "$progname: error at save table\n";

