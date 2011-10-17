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

if (defined($ARGV[0])&&("$ARGV[0]" eq "--help")){
  print <<EOF;
$progname, version %version%
This script:
- calculates the target function (ftar) by comparing the current and target density profile
- sorts ftars according to increasing order of magnitude in simplex table

Usage: $progname infile outfile param_N a_line_nr

NEEDS: name kBT

USES: readin_table csg_get_property saveto_table csg_get_property
EOF
  exit 0;
}

die "5 parameters are necessary\n" if ($#ARGV<4);

use CsgFunctions;
use SimplexFunctions;

my $infile="$ARGV[0]";
my $outfile="$ARGV[1]";
my $param_N="$ARGV[2]";
my $a_line_nr="$ARGV[3]";
my $prop_N="$ARGV[4]";

my $property="density";
my $name=csg_get_property("cg.non-bonded.name");
my $sim_prog=csg_get_property("cg.inverse.program");

# Resample tgt density
my $aim_dens_file="$name.dens.tgt";
my @args;
@args=("bash", "-c", "for_all non-bonded do_external resample_simplex $property");
system(@args);
undef(@args);

# Calculate new density
@args = ("bash", "-c", "for_all non-bonded do_external $property $sim_prog");
system(@args);
undef(@args);

my $cur_dens_file="$name.dens.new";

my @r_aim;
my @dens_aim;
my @flags_aim;
(readin_table($aim_dens_file,@r_aim,@dens_aim,@flags_aim)) || die "$progname: error at readin_table\n";

my @r_cur;
my @dens_cur;
my @flags_cur;
(readin_table($cur_dens_file,@r_cur,@dens_cur,@flags_cur)) || die "$progname: error at readin_table\n";

my @ftar_cur;
my @flag_cur;

my $ndim=$param_N+1;

# Read in temporary simplex table
my (%hash)=readin_simplex_table($infile,$ndim) or die "$progname: error at readin_simplex_table\n";

# Define table columns
@ftar_cur=@{$hash{p_0}};
@flag_cur=@{$hash{"p_$ndim"}};

# Should never happen due to resampling, but better check
die "Different grids \n" if (($r_aim[1]-$r_aim[0])!=($r_cur[1]-$r_cur[0]));
die "Different start point \n" if (($r_aim[0]-$r_cur[0]) > 0.0);

# ------------------- DEFINE TARGET FUNCTION HERE ------------------
# Calculate ftar
my @ddens=@_;
my $ftar=0;
my $ftar_aim=0;
my $dr=csg_get_interaction_property("inverse.simplex.density.step");
my $min=csg_get_interaction_property("inverse.simplex.density.min");
my $max=csg_get_interaction_property("inverse.simplex.density.max");

# Get absolute value
for(my $i=1;$i<($max-$min)/$dr;$i++) {
       $ddens[$i]=abs($dens_cur[$i]-$dens_aim[$i]);
       $ftar+=$ddens[$i];
       $ftar_aim+=$dens_aim[$i];
}

my $mdim;
if ($prop_N == 1) {
  $ftar_cur[$a_line_nr]=$ftar;
  $flag_cur[$a_line_nr]="complete";
  $mdim=$#ftar_cur+1;
}
else {
  $ftar_cur[0]=$ftar;
  for (my $j=1;$j<=$param_N;$j++){
    ${$hash{"p_$j"}}[0]=${$hash{"p_$j"}}[$a_line_nr];
  }
  $flag_cur[0]="complete";
  $mdim=1;
}

# Save to new simplex table
saveto_simplex_table($outfile,$mdim,$param_N,@ftar_cur,%hash,@flag_cur) or die "$progname: error at saveto_simplex_table\n";
