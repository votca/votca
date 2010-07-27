package CsgFunctions;
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

require Exporter;

use vars qw(@ISA @EXPORT);
@ISA         = qw(Exporter);
@EXPORT      = qw(csg_function_help csg_get_property csg_get_interaction_property readin_table readin_data saveto_table saveto_table_err);

sub csg_function_help() {
  print <<EOF;
CsgFunctions, version %version%
Provides useful function for perl:
csg_get_property($):             get a value from xml file
csg_get_interaction_property($): get a interaction property from xmlfile
readin_table(\$\\@\\@\\@):           reads in csg table
saveto_table(\$\\@\\@\\@):           writes to a csg table
saveto_table_err(\$\\@\\@\\@) :      writes to csg table with errors

USES: \$CSGXMLFILE csg_property
NEEDS:
PROVIDES: csg_get_property csg_get_interaction_property readin_table saveto_table saveto_table_err
EOF
  exit 0;
}

sub csg_get_property($){
  ( my $xmlfile=$ENV{'CSGXMLFILE'} ) || die "csg_get_property: ENV{'CSGXMLFILE'} was undefined\n";
  defined($_[0]) || die "csg_get_property: Missig argument\n";
  open(CSG,"csg_property --file $xmlfile --path $_[0] --short --print . |") ||
    die "csg_get_property: Could not open pipe\n";
  defined(my $value=<CSG>) || die "csg_get_property: Could not get value $_[0]\n";
  close(CSG) || die "csg_get_property: error from csg_property\n";
  chomp($value);
  return undef if ($value =~ /^\s*$/);
  return $value;
}

sub csg_get_interaction_property($){
  ( my $bondname=$ENV{'bondname'} ) || die "bondname: ENV{'bondname'} was undefined\n";
  ( my $bondtype=$ENV{'bondtype'} ) || die "bondtype: ENV{'bondtype'} was undefined\n";
  ( my $xmlfile=$ENV{'CSGXMLFILE'} ) || die "csg_get_property: ENV{'CSGXMLFILE'} was undefined\n";
  defined($_[0]) || die "csg_get_interaction_property: Missig argument\n";
  open(CSG,"csg_property --file $xmlfile --short --path cg.$bondtype --filter \"name=$bondname\" --print $_[0] |") ||
    die "csg_get_interaction_property: Could not open pipe\n";
  defined(my $value=<CSG>) || die "csg_get_interaction_property: Could not get value $_[0]\n";
  close(CSG) || die "csg_get_interaction_property: error from csg_property\n";
  chomp($value);
  return undef if ($value =~ /^\s*$/);
  return $value;
}

sub readin_table($\@\@\@) {
  defined($_[3]) || die "readin_table: Missing argument\n";
  open(TAB,"$_[0]") || die "readin_table: could not open file $_[0]\n";
  my $line=0;
  while (<TAB>){
    $line++;
    # remove leading spacees for split
    $_ =~ s/^\s*//;
    next if /^[#@]/;
    next if /^\s*$/;
    my @parts=split(/\s+/);
    defined($parts[1]) || die "readin_table: Not enought columns in line $line in file $_[0]\n";
    die "readin_table: Missing flag for r=$parts[0] in file $_[0]\n" unless (defined($parts[2]) && $parts[2] =~ /[iou]/);
    #very trick array dereference (@) of pointer to an array $_[.] stored in an array $_
    push(@{$_[1]},$parts[0]);
    push(@{$_[2]},$parts[1]);
    push(@{$_[3]},$parts[2]);
  }
  close(TAB) || die "readin_table: could not close file $_[0]\n";
  return $line;
}

sub readin_data($$\@\@) {
  defined($_[3]) || die "readin_data: Missing argument\n";
  open(TAB,"$_[0]") || die "readin_table: could not open file $_[0]\n";
  my $line=0;
  my $column=int($_[1]);
  while (<TAB>){
    $line++;
    # remove leading spacees for split
    $_ =~ s/^\s*//;
    next if /^[#@]/;
    next if /^\s*$/;
    my @parts=split(/\s+/);
    defined($parts[1]) || die "readin_table: Not enought columns in line $line in file $_[0]\n";
    die "readin_data: Can't read column $column\n" unless (defined($parts[$column]));
    #very trick array dereference (@) of pointer to an array $_[.] stored in an array $_
    push(@{$_[2]},$parts[0]);
    push(@{$_[3]},$parts[$column]);
  }
  close(TAB) || die "readin_table: could not close file $_[0]\n";
  return $line;
}

sub saveto_table($\@\@\@) {
  defined($_[3]) || die "saveto_table: Missing argument\n";
  open(OUTFILE,"> $_[0]") or die "saveto_table: could not open $_[0] \n";
  for(my $i=0;$i<=$#{$_[1]};$i++){
    print OUTFILE "${$_[1]}[$i] ${$_[2]}[$i] ${$_[3]}[$i]\n";
  }
  close(OUTFILE) or die "Error at closing $_[0]\n";
  return 1;
}

sub saveto_table_err($\@\@\@\@) {
  defined($_[3]) || die "saveto_table: Missing argument\n";
  open(OUTFILE,"> $_[0]") or die "saveto_table: could not open $_[0] \n";
  for(my $i=0;$i<=$#{$_[1]};$i++){
    print OUTFILE "${$_[1]}[$i] ${$_[2]}[$i] ${$_[3]}[$i] ${$_[4]}[$i]\n";
  }
  close(OUTFILE) or die "Error at closing $_[0]\n";
  return 1;
}

#important
1;
