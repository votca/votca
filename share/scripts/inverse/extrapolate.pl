#! /usr/bin/perl -w
#
# extrapolation methods:
#
# linear:
#   y = ax + b;  b = - m*xp + yp; a = m
# quadratic:
#   y = a*(x-b)^2; b = (x0 - 2y0/m); a = m^2/(4*y0)
# exponential:
#   y = a*exp(b*x); a = y0*exp(-m*x0/y0); b = m/y0;
#
#
#

sub extrapolate_linear($$$$) {
  my $x0 = $_[0];
  my $y0 = $_[1]; 
  my $m =  $_[2];
  my $x =  $_[3];

  return $m*($x - $x0) + $y0;
}

sub extrapolate_quad($$$$) {
  my $x0 = $_[0];
  my $y0 = $_[1]; 
  my $m =  $_[2];
  my $x =  $_[3];

  my $a = ($m**2)/(4*$y0);
  my $b = $x0 - 2*$y0/$m;
  #my $b = $x0 + 2*$y0/$m;
  #my $a = $m/(2*($x0-$b));

  return $a*($x-$b)**2;
}
sub extrapolate_exp($$$$) {
  my $x0 = $_[0];
  my $y0 = $_[1]; 
  my $m =  $_[2];
  my $x =  $_[3];

  my $a = $y0*exp(-$m*$x0 / $y0);
  my $b = $m/$y0;

  return $a*exp($b*$x);
}

use strict;

$_=$0;
s#^.*/##;
my $progname=$_;
my $usage="Usage: $progname [OPTIONS] <in> <out>";

# include functions
(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";


my $avgpoints = 3;
my $function="quadratic";
# read program arguments

while ((defined ($ARGV[0])) and ($ARGV[0] =~ /^\-/))
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
    if($ARGV[0] eq "--avgpoints") {
        $avgpoints = $ARGV[1];
        shift(@ARGV);
        shift(@ARGV);
    }
    elsif($ARGV[0] eq "--function") {
        $function = $ARGV[1];
        shift(@ARGV);
        shift(@ARGV);
    }
    elsif (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help"))
	{
		print <<END;
This script calculates the integral of a table
$usage
OPTIONS:
--avgpoints             average over so many points to extrapolate, standard is 3
--function            linear, quadratic or exp, standard is quadratic
-h, --help            Show this help message

END
		exit;
	}
	else
	{
		die "Unknow option '".$ARGV[0]."' !\n";
	}
}

#Print usage
die "no files given\n$usage\n" unless $#ARGV > 0;

my $infile="$ARGV[0]";
my @r;
my @val;
my @flag;
(readin_table($infile,\@r,\@val,\@flag)) || die "$progname: error at readin_table\n";

my $outfile="$ARGV[1]";

# find beginning
my $first;
for ($first=1;$first<=$#r;$first++) {
   last if($flag[$first] eq "i");
}

# find end
my $last;
for ($last=$#r;$last>0;$last--) {
   last if($flag[$last] eq "i");
}

# grad  of beginning
my $grad_beg;
$grad_beg = ($val[$first + $avgpoints] - $val[$first])/($r[$first + $avgpoints] - $r[$first]);

# grad  of beginning
my $grad_end;
$grad_end = ($val[$last] - $val[$last - $avgpoints])/($r[$last] - $r[$last-$avgpoints]);

# now extrapolate beginning
for(my $i=$first-1; $i >= 0; $i--) {
    $val[$i] = extrapolate_linear($r[$first], $val[$first], $grad_beg, $r[$i]);
}

# now extrapolate ends
for(my $i=$last+1; $i <= $#r; $i++) {
    $val[$i] = extrapolate_linear($r[$last], $val[$last], $grad_end, $r[$i]);
}

saveto_table($outfile,\@r,\@val,\@flag) || die "$progname: error at save table\n";
