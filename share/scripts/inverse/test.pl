#! /usr/bin/perl -w

use strict;

( my $progname = $0 ) =~ s#^.*/##;
(my $function_file=`$ENV{SOURCE_WRAPPER} functions perl`) || die "$progname: $ENV{SOURCE_WRAPPER} function perl failed\n";
chomp($function_file);
(do "$function_file") || die "$progname: source $function_file failed\n";

print csg_get("name");
