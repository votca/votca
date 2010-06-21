package SimplexFunctions;
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
use CsgFunctions;

require Exporter;

use vars qw(@ISA @EXPORT);
@ISA         = qw(Exporter);
@EXPORT      = qw(readin_init_simplex_table readin_simplex_table saveto_simplex_table sort_ftar calc_psum calc_ptry);

# Subroutine to read in simplex table
sub readin_init_simplex_table($$) {
  defined($_[1]) || die "readin_simplex_table: Missing file\n";
  my $infile=$_[0];
  my $param_N=$_[1];
  my %hash=();
  open(TAB,"$infile") || die "could not open file $_[0]\n";
  my $line=0;
  while (<TAB>) {
    $line++;
    # remove leading spaces for split
    $_ =~ s/^\s*//;
    next if /^[#@]/;
    next if /^\s*$/;
    my @values=split(/\s+/);
    defined($values[1]) || die "readin_table: Not enough columns in line $line in file $_[0]\n";
    foreach (1..$param_N) {
      push @{$hash{"p_$_"}}, $values[$_-1];
    }
  }
close(TAB) || die "could not close file $_[0]\n";
return %hash;
}

# Subroutine to read in simplex table
sub readin_simplex_table($$) {
  defined($_[1]) || die "readin_simplex_table: Missing file\n";
  my $infile=$_[0];
  my $ndim=$_[1];
  my %hash=();
  open(TAB,"$infile") || die "could not open file $_[0]\n";
  my $line=0;
  while (<TAB>) {
    $line++;
    # remove leading spaces for split
    $_ =~ s/^\s*//;
    next if /^[#@]/;
    next if /^\s*$/;
    my @values=split(/\s+/);
    defined($values[1]) || die "readin_table: Not enough columns in line $line in file $_[0]\n";
    foreach (0..$ndim) {
      push @{$hash{"p_$_"}}, $values[$_];
    }
  }
close(TAB) || die "could not close file $_[0]\n";
return %hash;
}

# Subroutine to save to simplex table
sub saveto_simplex_table($$$\@\%\@) {
  defined($_[5]) || die "saveto_table: Missing argument\n";
  my $outfile=$_[0];
  my $mdim=$_[1];
  my $param_N=$_[2];
  my @ftar=@{$_[3]};
  my (%hash)=%{$_[4]};
  my @flag=@{$_[5]};
  open(OUTFILE,"> $outfile") or die "could not open file $_[0]\n";
  for(my $i=0;$i<$mdim;$i++){
    print OUTFILE "$ftar[$i] ";
      for(my $j=1;$j<=$param_N;$j++){
        my @tmp=@{$hash{"p_$j"}};
        print OUTFILE "$tmp[$i] ";
      }
      print OUTFILE "$flag[$i]\n";
  }
  close(OUTFILE) or die "Error at closing $_[0]\n";
  return 1;
}

# Subroutine to sort columns of ftar
sub sort_ftar ($\@\@) {
  defined($_[2]) || die "get_psum: Missing argument: p[m][n]\n";
  my $param_N=$_[0];  
  my @ftar=@{$_[1]};
  my @p=@{$_[2]};  
  my $ndim=$param_N+1;
  my @i_sort;
  my @y_asc;
  my @p_asc;  
  foreach (0 .. $param_N) {
    $i_sort[$_]=$_;
  }
  @i_sort=(sort{$ftar[$a] <=> $ftar[$b]} @i_sort);
  for (my $i=0;$i<$ndim;$i++) {
    $y_asc[$i]=$ftar[$i_sort[$i]];
    for (my $j=0;$j<$param_N;$j++){
      $p_asc[$i][$j]=$p[$i_sort[$i]][$j];
    }
  }
  return (\@y_asc,\@p_asc);
}

# Subroutine to get sum of columns of p
sub calc_psum(\@$$) {
  defined($_[2]) || die "get_psum: Missing argument: p[m][n]\n";
  my @p=@{$_[0]};
  my $param_N=$_[1];
  my $ndim=$_[2];
  my $sum;
  my @psum;
  for(my $j=0; $j<$param_N;$j++) {
     for($sum=0, my $i=0;$i<$ndim;$i++) {
         $sum+=$p[$i][$j];
     }
     $psum[$j]=$sum;
     }
  return @psum;
}

# Subroutine calc_ptry
sub calc_ptry ($$$\@\@) {
defined($_[4]) || die "calc_ptry: Missing argument\n";
my $ndim="$_[0]";
my $ihi="$_[1]";
my $fac="$_[2]";
my @p=@{$_[3]};
my @psum=@{$_[4]};

my @ptry;
my $fac1=(1-$fac)/$ndim;
my $fac2=$fac1-$fac;

for (my $j=0;$j<$ndim;$j++) {
    $ptry[$j]=$psum[$j]*$fac1-$p[$ihi][$j]*$fac2;
    }
    return @ptry;
}

#important
1;

