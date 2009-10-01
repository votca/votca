#! /usr/bin/perl -w

use strict;
( my $progname = $0 ) =~ s#^.*/##;

sub find_from_table($);
sub find_in_dir($$);
sub search_and_exit;
sub show_table($);

if (defined($ARGV[0])&&( "$ARGV[0]" eq "--help" )){
  print <<END;
This script find a script from two keywords.
Usage: $progname word1 word2
First we check user table and search in CSGSCRIPTDIR and than in CSGINVERSE
Then we check csg table and search in CSGINVERSE

Hacky usage: $progname --direct scriptname
then we try to find scriptname in global,pwd,CSGSCRIPTDIR,CSGINVERSE

USES: \$CSGINVERSE \$CSGSCRIPTDIR 
NEEDS:
END
  exit 0;
}

my $csg_table="csg_table";
my $user_table="csg_table";

(my $csgshare=$ENV{'CSGINVERSE'}) || die "CSGINVERSE not defined\n";
$csg_table="$csgshare/$csg_table";

my $csgscriptdir=$ENV{'CSGSCRIPTDIR'};
if ($csgscriptdir) {
  $user_table="$csgscriptdir/$user_table";
} else {
  $user_table=undef;
  $csgscriptdir=undef;
}

my $scriptname=undef;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--status" )){
  print "csg table status\n";
  show_table($csg_table);
  print "Check sums\n";
  system('md5sum $CSGINVERSE/MD5SUM');
  system('cd $CSGINVERSE; md5sum -c $CSGINVERSE/MD5SUM || echo WARNING: You have modified csg scripts, better copy them and to user scripts dir');
  if (defined($user_table)&&( -r "$user_table")) {
    print "user table status\n";
    show_table($user_table);
  } else {
    print "No user table\n";
  }
  exit 0;
}

###################MAIN PROGRAMM#######################

( $#ARGV < 1 ) && die "$progname needs two arguments\n";

#--direct options
if ( "$ARGV[0]" eq "--direct" ) {
  search_and_exit($ARGV[1]," ",$ENV{'PWD'},$csgscriptdir,$csgshare);
  die "Could find script $ARGV[1]\n";
}

#find script in user_table
if (defined($user_table)&&( -r "$user_table" )&&($scriptname=find_from_table($user_table))) {
  #first script dir, then csgshare
  search_and_exit($scriptname,$csgscriptdir,$csgshare);
  die "Could not find user script '$scriptname' in any dir, check for typos\n";
}

#search in csg_table
( $scriptname=find_from_table($csg_table) ) || die "Could not find script matching '$ARGV[0] $ARGV[1]'\n";
search_and_exit($scriptname,$csgshare) || die "Could not find script '$scriptname' in any dir\n";

######################FUNCTIONS##########################
sub find_from_table($) {
  my $found=undef;
  open(FILE,"<$_[0]") || die "Could not read $_[0]\n";
  while (<FILE>) {
    next if /^#/;
    next if /^\s*$/;
    $found="$1" if /^\s*$ARGV[0]\s+$ARGV[1]\s+(.*?)$/;
  }
  return $found;
}

sub find_in_dir($$) {
  my $args="";
  ( my $script=$_[0] ) || die "find_in_dir: first argument missing\n";
  ( my $dir=$_[1] ) || die "find_in_dir: second argument missing\n";
  #remove script arguments
  $script="$1" if  ($_[0] =~ s/^(\S+).*?$/$1/);
  $args="$1" if ($_[0] =~ /^\S+\s+(.*?)$/); 
  #empty means, add no dir
  $script="$dir/$script"  unless ($dir =~ /^\s*$/ );
  if ( -f "$script" ) {
    $script.=" $args" unless ($args =~ /^\s*$/ );
    return "$script";
  }
  return undef;
}

sub search_and_exit {
  ( my $scriptname = shift ) || die "search_and_exit: first argument missing\n";
  foreach my $dir (@_) {
    next unless defined($dir);
    if ($_=find_in_dir($scriptname,$dir)) {
      print "$_\n";
      exit 0;
    }
  }
  return undef;
}

sub show_table($) {
  open(FILE,"<$_[0]") || die "Could not read $_[0]\n";
  while (<FILE>) {
    next if /^#/;
    next if /^\s*$/;
    print;
  }
  return 1;
}
