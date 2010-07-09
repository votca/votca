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
This script reads in the current and target PMFs, as well as 
the time-averaged particle distributions binned along a coordinate.
It outputs a file with the error made on each bead type.

Usage: $progname current_pmf target_pmf part_distribution num_bead_types outfile

NEEDS: cg.inverse.espresso.range_min cg.inverse.espresso.range_max cg.inverse.espresso.range_bins

USES: csg_get_property saveto_table readin_table readin_data
EOF
  exit 0;
}

die "5 parameters are nessary\n" if ($#ARGV<4);

use CsgFunctions;

my $cur_pmf    ="$ARGV[0]";
my $tgt_pmf    ="$ARGV[1]";
my $part_dist  ="$ARGV[2]";
my $num_p_types="$ARGV[3]";
my $outfile    ="$ARGV[4]";

die "Not enough bead types to analyze\n" if ($num_p_types < 1);

my $range_min=csg_get_property("cg.inverse.espresso.min");
my $range_max=csg_get_property("cg.inverse.espresso.max");
my $range_binsize=csg_get_property("cg.inverse.espresso.binsize");
my $pmf_ref_min=csg_get_property("cg.inverse.espresso.pmf_ref_min");
my $pmf_ref_max=csg_get_property("cg.inverse.espresso.pmf_ref_max");

# Current PMF
my @r_cur;
my @pmf_cur;
my @flag_cur;
(readin_table($cur_pmf,@r_cur,@pmf_cur,@flag_cur)) || die "$progname: error at readin_table\n";

# Target PMF
my @r_tgt;
my @pmf_tgt;
my @flag_tgt;
(readin_table($tgt_pmf,@r_tgt,@pmf_tgt,@flag_tgt)) || die "$progname: error at readin_table\n";

# Shift current and target PMFs. Reference interval used to set the zero.
# First calculate average value of reference interval in both curves.
my $avg_ref_cur=0.;
my $count_ref_cur=0;
my $avg_ref_tgt=0.;
my $count_ref_tgt=0;
for(my $i=0;$i<=$#r_cur;$i++){
    if ($r_cur[$i] < $pmf_ref_min || $r_cur[$i] > $pmf_ref_max) { 
    } else {
	# Current PMF
	$avg_ref_cur+=$pmf_cur[$i];
	$count_ref_cur++;
	# Target PMF
	$avg_ref_tgt+=$pmf_tgt[$i];
	$count_ref_tgt++;
    }
}
$avg_ref_cur/=$count_ref_cur;
$avg_ref_tgt/=$count_ref_tgt;
# Error on PMF
my @err_pmf;
for(my $i=0;$i<=$#r_cur;$i++){
    $pmf_cur[$i]-=$avg_ref_cur;
    $pmf_tgt[$i]-=$avg_ref_tgt;
    push(@err_pmf, abs($pmf_tgt[$i]/$pmf_cur[$i]));
}

my @update_factor;
for(my $t=1;$t<=$num_p_types;$t++){
    my @r_part;
    my @dist_part;
    (readin_data($part_dist, $t, @r_part, @dist_part)) || die "$progname: error at readin_data\n";

    # Total number of particles of type $t
    my $n_part_type=0;
    ($n_part_type+=$_) for @dist_part;

    # Cumulative error on PMF
    my $cumulative_err=0.;


    printf("***************\ntype #$t - $n_part_type beads\n");
    printf("r part_dist | cur_pmf tgt_pmf | error\n");
    for(my $i=0;$i<=$#r_part;$i++){
	$cumulative_err+= ($dist_part[$i]/$n_part_type) * $err_pmf[$i];
	printf("%f %f | %f %f | %f\n", $r_part[$i], $dist_part[$i], $pmf_cur[$i], $pmf_tgt[$i], $err_pmf[$i]);
    }
    printf("cumul error: %f\n",$cumulative_err);
    push(@update_factor,$cumulative_err);
}


# Now output update factor for all bead types
open(OUTFILE,"> $outfile") or die "saveto_table: could not open $outfile\n";
for(my $t=0;$t<$num_p_types;$t++){
  printf(OUTFILE "%f", $update_factor[$t]);
  if ($t<$num_p_types-1) {
      printf(OUTFILE " ");
  } else {
      printf(OUTFILE "\n");
  }
}
close(OUTFILE) or die "Error at closing $outfile\n";

