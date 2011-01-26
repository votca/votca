#! /usr/bin/perl -w
#
# Copyright 2009 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

use strict;
( my $progname = $0 ) =~ s#^.*/##;

sub find_from_table($);
sub find_in_dir($$);
sub search_and_exit;
sub show_table($);

if (defined($ARGV[0])&&( "$ARGV[0]" eq "--help" )){
  print <<END;
$progname, version %version%
This script find a script from two keywords.
First we check user table and search in CSGSCRIPTDIR and than in CSGINVERSE
Then we check csg table and search in CSGINVERSE

Hacky usage: $progname --direct scriptname
then we try to find scriptname in global,pwd,CSGSCRIPTDIR,CSGINVERSE

Usage: $progname word1 word2

USES: \$CSGINVERSE \$CSGSCRIPTDIR

END
  exit 0;
}

my $csg_table="csg_table";
my $user_table="csg_table";

(my $csgshare=$ENV{'CSGINVERSE'}) or die "CSGINVERSE not defined\n";
$csg_table="$csgshare/$csg_table";
( -r "$csg_table" ) or die "Could read $csg_table\n";

my $csgscriptdir=$ENV{'CSGSCRIPTDIR'};
if ($csgscriptdir) {
  $user_table="$csgscriptdir/$user_table";
  unless ( -r "$user_table" ) {
    $user_table=undef;
    $csgscriptdir=undef;
  }
} else {
  $user_table=undef;
  $csgscriptdir=undef;
}

my $scriptname=undef;

if (defined($ARGV[0])&&("$ARGV[0]" eq "--status" )){
  print "csg table status\n";
  show_table($csg_table);
  if (defined($user_table)) {
    print "user table status (userscriptdir: $csgscriptdir)\n";
    show_table($user_table);
  } else {
    print "No user table\n";
  }
  exit 0;
}

if (defined($ARGV[0])&&("$ARGV[0]" eq "--check" )){
  print "Check sums\n";
  ( -r "$csgshare/MD5SUM" ) or die "Could not read checksum file $csgshare/MD5SUM\n";
  system('md5sum $CSGINVERSE/MD5SUM');
  system('cd $CSGINVERSE;md5sum -c MD5SUM || echo WARNING: You have modified csg scripts, better copy them to your user scripts dir');
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
if (defined($user_table)&&($scriptname=find_from_table($user_table))) {
  #first script dir, then csgshare
  search_and_exit($scriptname,$csgscriptdir,$csgshare);
  die "Could not find script '$ARGV[0] $ARGV[1]' ($scriptname) in user and csg dir, check for typos\n";
}

#search in csg_table
( $scriptname=find_from_table($csg_table) ) or die "Could not find script matching '$ARGV[0] $ARGV[1]'\n";
search_and_exit($scriptname,$csgshare) or die "Could not find script '$ARGV[0] $ARGV[1]' ($scriptname) in csg dir\n";

######################FUNCTIONS##########################
sub find_from_table($) {
  my $found=undef;
  open(FILE,"<$_[0]") or die "Could not read $_[0]\n";
  while (<FILE>) {
    next if /^#/;
    next if /^\s*$/;
    $found="$1" if /^\s*$ARGV[0]\s+$ARGV[1]\s+(.*?)$/;
  }
  return $found;
}

sub find_in_dir($$) {
  my $args="";
  ( my $script=$_[0] ) or die "find_in_dir: first argument missing\n";
  ( my $dir=$_[1] ) or die "find_in_dir: second argument missing\n";
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
  ( my $scriptname = shift ) or die "search_and_exit: first argument missing\n";
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
  open(FILE,"<$_[0]") or die "Could not read $_[0]\n";
  while (<FILE>) {
    next if /^#/;
    next if /^\s*$/;
    print;
  }
  return 1;
}
