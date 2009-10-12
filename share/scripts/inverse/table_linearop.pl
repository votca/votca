#! /usr/bin/perl -w
use strict;

$_=$0;
s#^.*/##;
my $progname=$_;
my $usage="Usage: $progname [OPTIONS] <in> <out> <a> <b>";

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
Performs a linear operaton on the y values:
y_new = a*y_old + b
$usage
OPTIONS:
-h, --help            Show this help message
--withflag            only change entries with specific flag in src

Examples:  $progname tmp.dpot.cur tmp.dpot.new 1.0 0.0

USES: \$SOURCE_WRAPPER readin_table saveto_table
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
die "missing parameters\n$usage\n" unless $#ARGV >= 3;

my $a = $ARGV[1];
my $b = $ARGV[2]; 

# include perl functions
(my $function_file=`$ENV{'SOURCE_WRAPPER'} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";
########

my $file="$ARGV[0]";

print "table $file : y' = $a*y + $b\n";

my @r;
my @val;
my @flag;
(readin_table($file,\@r,\@val,\@flag)) || die "$progname: error at readin_table\n";

for(my $i=0; $i<=$#r; $i++) {
  # skip if flag does not match
  if($withflag) {
    if(!($flag[$i] =~ m/[$withflag]/)) {
      next;
    }
  }
  $val[$i] = $a*$val[$i] + $b;
}

saveto_table($file,\@r,\@val,\@flag) || die "$progname: error at save table\n";

