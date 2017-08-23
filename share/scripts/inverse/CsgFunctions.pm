package CsgFunctions;
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

require Exporter;

use vars qw(@ISA @EXPORT);
@ISA         = qw(Exporter);
@EXPORT      = qw(csg_function_help readin_table readin_data saveto_table saveto_table_err readin_table_err);

sub csg_function_help() {
  print <<EOF;
CsgFunctions, version %version%
Provides useful function for VOTCA's iterative framework written in perl:
readin_table(\$\\@\\@\\@;\\\$):           reads in csg table
readin_table_err(\$\\@\\@\\@;\\\$):           reads in csg table with errors
saveto_table(\$\\@\\@\\@;\$):           writes to a csg table
saveto_table_err(\$\\@\\@\\@;\$) :      writes to csg table with errors
EOF
  exit 0;
}

sub readin_table($\@\@\@;\$) {
  defined($_[3]) || die "readin_table: Missing argument\n";
  open(TAB,"$_[0]") || die "readin_table: could not open file $_[0]\n";
  my $sloppy= $ENV{'VOTCA_TABLES_WITHOUT_FLAG'};
  $sloppy="no" unless defined($sloppy);
  my $line=0;
  while (<TAB>){
    $line++;
    ${$_[4]}.=$_ if (defined($_[4]) and (/^[#@]/));
    next if /^[#@]/;
    # remove leading spacees for split
    $_ =~ s/^\s*//;
    next if /^\s*$/;
    my @parts=split(/\s+/);
    if ( $sloppy eq "yes" ) {
      defined($parts[1]) || die "readin_table: Not enought columns in line $line in file $_[0]\n";
      $parts[$#parts+1] = "i";
    } else {
      defined($parts[2]) || die "readin_table: Not enought columns in line $line in file $_[0], if you don't have flags in your table add --sloppy-tables option to csg_call\n";
      ($parts[$#parts] =~ /[iou]/) || die "readin_table: Wrong flag($parts[$#parts]) for r=$parts[0] in file $_[0], if you don't have flags in your table add --sloppy-tables option to csg_call\n";
    }
    #very trick array dereference (@) of pointer to an array $_[.] stored in an array $_
    push(@{$_[1]},$parts[0]);
    push(@{$_[2]},$parts[1]);
    push(@{$_[3]},$parts[$#parts]);
  }
  close(TAB) || die "readin_table: could not close file $_[0]\n";
  die "readin_table: 0 lines were read from $_[0]\n" if ($line==0);
  return $line;
}

sub readin_table_err($\@\@\@\@;\$) {
  defined($_[4]) || die "readin_table_err: Missing argument\n";
  open(TAB,"$_[0]") || die "readin_table_err: could not open file $_[0]\n";
  my $sloppy= $ENV{'VOTCA_TABLES_WITHOUT_FLAG'};
  $sloppy="no" unless defined($sloppy);
  my $line=0;
  while (<TAB>){
    $line++;
    ${$_[5]}.=$_ if (defined($_[4]) and (/^[#@]/));
    # remove leading spacees for split
    $_ =~ s/^\s*//;
    next if /^[#@]/;
    next if /^\s*$/;
    my @parts=split(/\s+/);
    if ( $sloppy eq "yes" ) {
      defined($parts[2]) || die "readin_table_err: Not enought columns in line $line in file $_[0]\n";
      $parts[$#parts+1] = "i";
    }else{
      defined($parts[3]) || die "readin_table_err: Not enought columns in line $line in file $_[0], if you don't have flags in your table add --sloppy-tables option to csg_call\n";
      ($parts[$#parts] =~ /[iou]/) || die "readin_table_err: Wrong flag($parts[$#parts]) for r=$parts[0] in file $_[0], if you don't have flags in your table add --sloppy-tables option to csg_call\n";
    }
    #very trick array dereference (@) of pointer to an array $_[.] stored in an array $_
    push(@{$_[1]},$parts[0]);
    push(@{$_[2]},$parts[1]);
    push(@{$_[3]},$parts[2]);
    push(@{$_[4]},$parts[$#parts]);
  }
  close(TAB) || die "readin_table_err: could not close file $_[0]\n";
  die "readin_table_err: 0 lines were read from $_[0]\n" if ($line==0);
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

sub saveto_table($\@\@\@;$) {
  defined($_[3]) || die "saveto_table: Missing argument\n";
  open(OUTFILE,"> $_[0]") or die "saveto_table: could not open $_[0] \n";
  print OUTFILE "$_[4]" if (defined $_[4]);
  for(my $i=0;$i<=$#{$_[1]};$i++){
    print OUTFILE "${$_[1]}[$i] ${$_[2]}[$i] ${$_[3]}[$i]\n";
  }
  close(OUTFILE) or die "Error at closing $_[0]\n";
  return 1;
}

sub saveto_table_err($\@\@\@\@;$) {
  defined($_[3]) || die "saveto_table: Missing argument\n";
  open(OUTFILE,"> $_[0]") or die "saveto_table: could not open $_[0] \n";
  print OUTFILE "$_[5]" if (defined $_[5]);
  for(my $i=0;$i<=$#{$_[1]};$i++){
    print OUTFILE "${$_[1]}[$i] ${$_[2]}[$i] ${$_[3]}[$i] ${$_[4]}[$i]\n";
  }
  close(OUTFILE) or die "Error at closing $_[0]\n";
  return 1;
}

#important
1;
