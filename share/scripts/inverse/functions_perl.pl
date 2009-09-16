#! /usr/bin/perl -w

use strict;

if ("$ARGV[0]" eq "--help"){
  print <<EOF;
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
    unless (defined ($parts[2])){
      $parts[2]="?";
    } 
    defined($parts[1]) || die "readin_table: Not enought columns in line $line\n";
    #very trick array dereference (@) of pointer to an array $_[.] stored in an array $_
    push(@{$_[1]},$parts[0]);
    push(@{$_[2]},$parts[1]);
    push(@{$_[3]},$parts[2]);
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
return 1;
