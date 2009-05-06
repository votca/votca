#! /usr/bin/perl -w

use strict;

sub get_sim_property($){
  ($_[0]) || die "get_sim_property: Missig argument\n";
  open(CSG,"csg_property --file $ENV{'CSGXMLFILE'} --path cg.inverse.$_[0] --short --print . |") || 
    die "get_sim_property: Could not open pipe\n";
  defined(my $value=<CSG>) || die "get_sim_property: Could not get value $_[0]\n";
  close(CSG) || die "get_sim_property: error from csg_property\n";
  chomp($value);
  return $value;
}

sub csg_get($){
  ($_[0]) || die "csg_get: Missig argument\n";
  open(CSG,"$ENV{'csg_get'} $_[0] |") || die "csg_get: Could not open pipe\n";
  defined(my $value=<CSG>) || die "csg_get: Could not get value $_[0]\n";
  close(CSG) || die "csg_get: error from csg_property\n";
  chomp($value);
  return $value;
}

sub readin_table($\@\@\@) {
  ($_[3]) || die "readin_table: Missing argument\n";
  open(TAB,"$_[0]") || die "readin_table: could not open file $_[0]\n";
  my $line=0;
  while (<TAB>){
    $line++;
    next if /^[#@]/;
    next if /^\s*$/;
    my @parts=split(/\s+/);
    ($parts[2]) || die "readin_table: Not enought columns in line $line\n";
    #very trick array dereference (@) of pointer to an array $_[.] stored in an array $_
    push(@{$_[1]},$parts[0]);
    push(@{$_[2]},$parts[1]);
    push(@{$_[3]},$parts[2]);
  }
  close(TAB) || die "readin_table: could not close file $_[0]\n";
  return $line;
}

sub saveto_table($\@\@\@) {
  ($_[3]) || die "saveto_table: Missing argument\n";
  open(OUTFILE,"> $_[0]") or die "saveto_table: could not open $_[0] \n";
  for(my $i=0;$i<=$#{$_[1]};$i++){
    print OUTFILE "${$_[1]}[$i] ${$_[2]}[$i] ${$_[3]}[$i]\n";
  }
  close(OUTFILE) or die "Error at closing $_[0]\n";
  return 1;
}

#important
return 1;
