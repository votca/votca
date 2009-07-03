#! /usr/bin/perl -w
#
# (C) 2006-2008 Chr. Junghans
# junghans@mpip-mainz.mpg.de
#
#
#version 0.1  , 08.07.08 -- initial version
use strict;

$_=$0;
s#^.*/##;
my $progname=$_;
my $usage="Usage: $progname [OPTIONS] <in> <out>";

# include functions
(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";

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
	if (($ARGV[0] eq "-h") or ($ARGV[0] eq "--help"))
	{
		print <<END;
This script calculates the integral of a table
$usage
OPTIONS:
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
my @out;

my $min = 0;
$out[0] = 0;
for (my $i=1;$i<=$#r;$i++){
  $out[$i]=$out[$i-1] - 0.5*($val[$i] + $val[$i-1])*($r[$i] - $r[$i-1]); 
  $min = $out[$i] if $out[$i] < $min;
}

for (my $i=0;$i<=$#r;$i++){
  $out[$i]=$out[$i]-$min;
}

saveto_table($outfile,\@r,\@out,\@flag) || die "$progname: error at save table\n";
