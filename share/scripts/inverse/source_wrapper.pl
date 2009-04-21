#! /usr/bin/perl -w

use strict;
( my $progname = $0 ) =~ s#^.*/##;

sub find_from_table($);
sub find_in_dir($$);
sub search_and_exit;

if ( "$ARGV[0]" eq "--help" ){
  print <<END;
This script find a script from two keywords.
Usage: $progname word1 word2
First we check user table and search in CSGSCRIPTDIR and than in CSGSHARE
Then we check csg table and search in CSGSHARE

Hacky usage: $progname --direct scriptname
then we try to find scriptname in global,pwd,CSGSCRIPTDIR,CSGSHARE
END
  exit 0;
}

my $csg_table="csg_table";
my $user_table="csg_table";

(my $csgshare=$ENV{'CSGSHARE'}) || die "CSGSHARE not defined\n";
(my $csgscriptdir=$ENV{'CSGSCRIPTDIR'}) ||die "CSGSCRIPTDIR not defined\n";

$csg_table="$csgshare/$csg_table";
$user_table="$csgscriptdir/$user_table";

my $scriptname=undef;

###################MAIN PROGRAMM#######################

( $#ARGV < 1 ) && die "$progname needs two arguments\n";

#--direct options
if ( "$ARGV[0]" eq "--direct" ) {
  search_and_exit($ARGV[1]," ",$ENV{'PWD'},$csgscriptdir,$csgshare);
  die "Could find script $ARGV[1]\n";
}

#find script in user_table
if (( -r "$user_table" ) && ($scriptname=find_from_table($user_table))) {
  #first script dir, then csgshare
  search_and_exit($scriptname,$csgscriptdir,$csgshare);
}

#search in csg_table
( $scriptname=find_from_table($csg_table) ) || die "Could not find script matching '$ARGV[0] $ARGV[1]'\n";
search_and_exit($scriptname,$csgshare) || die "Could find script '$scriptname'\n";

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
    if ($_=find_in_dir($scriptname,$dir)) {
      print "$_\n";
      exit 0;
    }
  }
  return undef;
}
