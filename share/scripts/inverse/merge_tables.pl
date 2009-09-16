#! /usr/bin/perl -w
use strict;

$_=$0;
s#^.*/##;
my $progname=$_;
my $usage="Usage: $progname [OPTIONS] <source> <dest> <out>";

#Defaults
my $noflags='no';
my $novalues='no';
my $withflag=undef;

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
Merge two tables
$usage
OPTIONS:
-v, --version         Prints version
-h, --help            Show this help message
--withflag            only change entries with specific flag in src
--noflags             don't copy flags
--novalues            don't copy values

Examples:  $progname -q
           $progname

USES: \$SOURCE_WRAPPER readin_table saveto_table
NEEDS: 

END
		exit;
	}
    elsif ($ARGV[0] eq "--noflags")
    {
        $noflags = 'yes'
    }
    elsif ($ARGV[0] eq "--novalues")
    {
        $novalues = 'yes'
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
die "missing parameters\n$usage\n" unless $#ARGV > 1;

# include perl functions
(my $function_file=`$ENV{'SOURCE_WRAPPER'} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";
########

my $src="$ARGV[0]";
my $dst="$ARGV[1]";
my $out="$ARGV[2]";

print "tables $src $dst $out\n";

my @r_src;
my @val_src;
my @flag_src;
(readin_table($src,\@r_src,\@val_src,\@flag_src)) || die "$progname: error at readin_table\n";

my @r_dst;
my @val_dst;
my @flag_dst;
(readin_table($dst,\@r_dst,\@val_dst,\@flag_dst)) || die "$progname: error at readin_table\n";

my $idst=0;

for(my $i=0; $i<=$#r_src; $i++) {
  # skip if flag does not match
  if($withflag) {
    if(!($flag_src[$i] =~ m/[$withflag]/)) {
      next;
    }
  }
  
  # advance in dst till same r
  while($r_dst[$idst] < $r_src[$i] - 1e-15) {
    $idst++;
  }
  my $tmp= $r_src[$i]-$r_dst[$idst];

  die "error: grid mismatch" if(abs($r_dst[$idst] - $r_src[$i]) > 1e-15);

  if($novalues eq 'no') {
    $val_dst[$idst] = $val_src[$i];
  }
  if($noflags eq 'no') {
    $flag_dst[$idst] = $flag_src[$i];
  }
}

saveto_table($out,\@r_dst,\@val_dst,\@flag_dst) || die "$progname: error at save table\n";

