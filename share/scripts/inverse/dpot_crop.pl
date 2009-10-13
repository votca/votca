#! /usr/bin/perl -w
use strict;

$_=$0;
s#^.*/##;
my $progname=$_;
my $usage="Usage: $progname [OPTIONS] <file> <a> <b>";

#Defaults
my $withflag=undef;

while ((defined ($ARGV[0])) and ($ARGV[0] =~ /^-./))
{
        if (($ARGV[0] !~ /^--/) and (length($ARGV[0])>2)){
           $_=shift(@ARGV);
           #short opt having agruments examples fo
           if ( $_ =~ /^-[fo]/ ) {
              unshift(@ARGV,substr($_,0,2),substr($_,2));
           }
           else{
              unshift(@ARGV,substr($_,0,2),"-".substr($_,2));
           }
        }
	if (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help"))
	{
		print <<END;
crop the potential update at poorly sampled ends
$usage
OPTIONS:
-h, --help            Show this help message

Examples:  $progname tmp.dpot.cur tmp.dpot.new 

USES: readin_table saveto_table
NEEDS: 

END
		exit;
	}
    elsif ($ARGV[0] eq "--withflag")
    {
        shift(@ARGV);
        die "nothing given for --withflag" unless $#ARGV > -1;
        $withflag = $ARGV[0];
    }
	else
	{
		die "Unknow option '".$ARGV[0]."' !\n";
	}
    shift(@ARGV);
}

#Print usage
die "missing parameters\n$usage\n" unless $#ARGV >= 1;

use CsgFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";

my @r;
my @val;
my @flag;
(readin_table($infile,@r,@val,@flag)) || die "$progname: error at readin_table\n";

# find last u/o
my $i_first;

# TODO: look for at least 3 successive points with i
for($i_first=0; ($i_first<$#r) && ($flag[$i_first] =~ /[uo]/); $i_first++) {}

my $ncrop=0;

while($i_first + $ncrop<=$#r-3) {
  my $i = $i_first + $ncrop;
  my $delta_1 = $val[$i] -  $val[$i + 1];
  my $delta_2 = $val[$i + 1 ] -  $val[$i + 2];

  # do both deltas have the same sign?
  if($delta_1 * $delta_2 > 0) {
    last;
  } elsif (abs($val[$i]) < 0.5 && abs($val[$i+1]) < 0.5) {
    last; 
  }
  $flag[$i]='o';
  $ncrop++;
  if($ncrop > 3) {
    print "error: need to crop more than 3 points in $infile. think about sampleing/grid interval.";
    exit 1;
  }
}

if($ncrop > 0) {
  print "warnng, I cropped $ncrop points at the beginning\n";
}

saveto_table($outfile,@r,@val,@flag) || die "$progname: error at save table\n";

