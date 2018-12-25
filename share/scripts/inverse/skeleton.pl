#! /usr/bin/perl -w
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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
if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
$progname, version %version%
Skeleton script

Usage: $progname argument
EOF
  exit 0;
}

die "1 parameter are necessary\n" if ($#ARGV<0);

use CsgFunctions;

my @r;
my @values;
my @flags;
readin_table($ARGV[0],@r,@values,@flags)) || die "$progname: error at readin_table'$ARGV[0]'\n"
