#! /usr/bin/perl -w

use strict;

sub get_sim_property{
  ($_[0]) || die "get_sim_property: Missig arrgument for get_sim_property";
  open(CSG,"csg_property --file $ENV{'CSGXMLFILE'} --path cg.inverse.$_[0] --short --print |") || 
    die "get_sim_property: Could not open pipe\n";
  defined(my $value=<CSG>) || die "get_sim_property: Could not get value $_[0]\n";
  close(CSG) || die "get_sim_property: error from csg_property\n";
  return $_;
}

sub csg_get{
  ($_[0]) || die "csg_get: Missig arrgument for get_sim_property";
  open(CSG,"$ENV{'csg_get'} $_[0] |") || die "csg_get: Could not open pipe\n";
  defined(my $value=<CSG>) || die "csg_get: Could not get value $_[0]\n";
  close(CSG) || die "csg_get: error from csg_property\n";
  return $value;
}

#important
return 1;
