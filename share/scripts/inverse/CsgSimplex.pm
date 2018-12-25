package CsgSimplex;
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
@EXPORT      = qw(csg_simplex_function_help readin_simplex_state saveto_simplex_state sort_simplex_table is_num replace_parameter_flag get_convergence_value remove_parameter_set calc_parameter_center linop_parameter);

sub csg_simplex_function_help() {
  print <<EOF;
CsgSimplex, version %version%
Provides useful simplex function for VOTCA's iterative framework in perl:
readin_simplex_state(\$\\\$\\@;\\\$):  reads in simplex state
saveto_simplex_state(\$\$\\@;\$):      writes simplex state
sort_simplex_table(\\@);               sort simplex parametertable
is_num(\$);
EOF
  exit 0;
}

sub readin_simplex_state($\$\@;\$\$) {
  defined($_[2]) || die "readin_simplex_state: Missing argument\n";
  open(TAB,"$_[0]") || die "readin_simplex_state: could not open file $_[0]\n";
  my $line=0;
  my $parameters=undef;
  my $state=undef;
  while (<TAB>){
    $line++;
    if (/^#State =\s*(\S*)/) {
      $state="$1";
      next;
    }
    ${$_[3]}.=$_ if (defined($_[3]) and (/^[#@]/));
    if (/^#Format\s*(.*)$/) {
      ${$_[4]}="$1" if defined($_[4]);
    }
    next if ($_ =~ /^[#@]/);
    $_ =~ s/^\s*//; # remove leading spacees for split
    next if /^\s*$/;
    my @parts=split(/\s+/);
    $parameters=$#parts unless ($parameters);
    die "readin_simplex_state: Number of parameters ($#parts) differ in line $line from previous lines ($parameters) of file $_[0]\n" if ($parameters != $#parts);
    push(@{$_[2]},\@parts);
  }
  close(TAB) || die "readin_simplex_state: could not close file $_[0]\n";
  die "readin_simplex_state: 0 lines were read from $_[0]\n" if ($line==0);
  die "readin_simplex_state: could not read state from  $_[0]\n" unless ($state);
  ${$_[1]}="$state";
  return $line;
}

sub saveto_simplex_state($$\@;$) {
  defined($_[2]) || die "saveto_simplex_state: Missing argument\n";
  open(OUTFILE,"> $_[0]") or die "saveto_simplex_state: could not open $_[0] \n";
  print OUTFILE "$_[3]" if (defined $_[3]);
  print OUTFILE "#State = $_[1]\n";
  my @simplex_table=@{$_[2]};
  my $parameters=undef;
  #remember 2d arrays is a list of lists
  for (my $i=0;$i<=$#simplex_table;$i++){
    $parameters=$#{$simplex_table[$i]} unless ($parameters);
    die "saveto_simplex_state: Number of parameters ($#{$simplex_table[$i]}) differ in set $i from previous lines ($parameters)"  if ($parameters != $#{$simplex_table[$i]});
    print OUTFILE "@{$simplex_table[$i]}\n";
  }
  close(OUTFILE) or die "Error at closing $_[0]\n";
  return 1;
}

sub sort_simplex_table(\@) {
  defined($_[0]) || die "sort_simplex_table: Missing argument\n";
  my @simplex_table=@{$_[0]};
  my $parameters=undef;
  my @index;
  #remember 2d arrays is a list of lists
  for (my $i=0;$i<=$#simplex_table;$i++){
    $parameters=$#{$simplex_table[$i]} unless ($parameters);
    die "sort_simplex_table: Number of parameters ($#{$simplex_table[$i]}) differ in set $i from previous lines ($parameters)"  if ($parameters != $#{$simplex_table[$i]});
    push(@index,$i);
    is_num("$simplex_table[$i][-2]") || die "sort_simplex_table: second last value of parameter set $i is not a number\n";
    die "sort_simplex_table: try set found!\n" if ($simplex_table[$i][-1] =~ /^try$/);
  }
  @index=sort { $simplex_table[$a][-2] <=> $simplex_table[$b][-2] } @index;
  my @sorted_table;
  for (my $i=0;$i<=$#index;$i++){
    $sorted_table[$i]=$simplex_table[$index[$i]];
  }
  @{$_[0]}=@sorted_table;
  return 1;
}

sub replace_parameter_flag(\@$$) {
  defined($_[2]) || die "replace_parameter_flag: Missing argument\n";
  my @simplex_table=@{$_[0]};
  for (my $i=0;$i<=$#simplex_table;$i++){
    $simplex_table[$i][$#{$simplex_table[$i]}] =~ s/^$_[1]$/$_[2]/;
  }
}


sub is_num($) {
  defined($_[0]) || die "is_num: Missing argument\n";
  use POSIX qw(strtod);
  my $str = shift;
  $str =~ s/^\s+//;
  $str =~ s/\s+$//;
  $! = 0;
  my($num, $unparsed) = strtod($str);
  if (($str eq '') || ($unparsed != 0) || $!) { return 0; }
  return 1;
}

# prototype is needed here as the function calls itself.
sub get_convergence_value(\@$);
sub get_convergence_value(\@$) {
  defined($_[1]) || die "get_convergence_value: Missing argument\n";
  my @simplex_table=@{$_[0]};
  if ($_[1] eq "lowest") {
    my $value=undef;
    for (my $i=0;$i<=$#simplex_table;$i++) {
      next if ($simplex_table[$i][-1] =~ /^try/);
      $value=$simplex_table[$i][-2] unless defined($value);
      $value=$simplex_table[$i][-2] if $value>$simplex_table[$i][-2]
    }
    return $value;
  } elsif ($_[1] eq "ihighest") {
    my $ivalue=undef;
    for (my $i=0;$i<=$#simplex_table;$i++) {
      next if ($simplex_table[$i][-1] =~ /^(try|tryold)$/);
      $ivalue=$i unless defined($ivalue);
      $ivalue=$i if $simplex_table[$ivalue][-2]<$simplex_table[$i][-2];
    }
    return $ivalue;
  } elsif ($_[1] eq "highest") {
    my $i=get_convergence_value(@simplex_table,"ihighest");
    return $simplex_table[$i][-2];
  } elsif ($_[1] eq "second") {
    my $ivalue=get_convergence_value(@simplex_table,"ihighest");
    # in case we do simplex on one parameter
    return $simplex_table[$ivalue][-2] if ($#simplex_table == 2);
    my $value=undef;
    for (my $i=0;$i<=$#simplex_table;$i++) {
      next if ($simplex_table[$i][-1] =~ /^(try|tryold)$/);
      next if ($i==$ivalue);
      $value=$simplex_table[$i][-2] unless defined($value);
      $value=$simplex_table[$i][-2] if $value<$simplex_table[$i][-2];
    }
    return $value;
  } elsif ($_[1] =~ /^(try|tryold)$/) {
    my $value=undef;
    for (my $i=0;$i<=$#simplex_table;$i++) {
      if ( $simplex_table[$i][-1] =~ /^$_[1]$/ ) {
        die "get_convergence_value: Found two $_[1] value in parameter set\n" if defined($value);
        $value=$simplex_table[$i][-2];
      }
    }
    die "get_convergence_value: Could not find any $_[1] value in paramter set\n" unless defined($value);
    return $value;
  } else { 
    die "get_convergence_value: I don't know how to get value '$_[1]'\n";
  }
}

sub remove_parameter_set(\@$) {
  defined($_[1]) || die "remove_parameter_set: Missing argument\n";
  my @simplex_table=@{$_[0]};
  my $value=undef;
  my @new_table;
  if ($_[1] =~ /^(try|tryold)$/) {
    for (my $i=0;$i<=$#simplex_table;$i++) {
      if ( $simplex_table[$i][-1] =~ /^$_[1]$/ ) {
        die "remove_parameter_set: Found two parameter set with flag '$_[1]'" if defined($value);
        $value=$i;
      } else {
        push(@new_table,$simplex_table[$i]);
      }
    }
  } elsif ($_[1] eq "highest") {
    $value=get_convergence_value(@simplex_table,"ihighest");
    for (my $i=0;$i<=$#simplex_table;$i++) {
      push(@new_table,$simplex_table[$i]) unless ($i == $value);
    }
  } else {
    die "remove_parameter_set: I don't know how to remove value '$_[1]'\n";
  }
  die "remove_parameter_set: Could not find a parameter set with flag '$_[1]'" unless defined($value);
  @{$_[0]}=@new_table;
  return @{$simplex_table[$value]};
}

sub calc_parameter_center(\@){
  defined($_[0]) || die "calc_parameter_center: Missing argument\n";
  my @simplex_table=@{$_[0]};
  my @center;
  sort_simplex_table(@simplex_table);
  for (my $j=0;$j<=$#{$simplex_table[0]};$j++) {
    if (is_num("$simplex_table[0][$j]")) {
      $center[$j]=0;
    } else {
      $center[$j]=$simplex_table[0][$j];
    }
  }
  #mind the $i<$#simplex_table to skip the highest value
  for (my $i=0;$i<$#simplex_table;$i++) {
    die "calc_parameter_center: number of parameters (".($#{$simplex_table[$i]}-1).") of parameter set #".($i+1)." differs from the number of non-try sets - 1 (=".($#simplex_table)."). Expected $#{$simplex_table[$i]} non-try sets.\n" if (($#simplex_table+1) != $#{$simplex_table[$i]});

    for (my $j=0;$j<=$#{$simplex_table[$i]};$j++) {
      if (is_num("$simplex_table[$i][$j]")) {
	$center[$j]+=$simplex_table[$i][$j]/$#simplex_table;
      }
    }
  }
  return @center;
}

sub linop_parameter(\@$\@\@) {
  defined($_[3]) || die "linop_parameter: Missing argument\n";
  my @vec1=@{$_[0]};
  my $scale=$_[1];
  my @vec2=@{$_[2]};
  die "linop_parameter: 1st ($#vec1) and 2nd vector ($#vec2) have different length\n" unless ($#vec1 = $#vec2);
  my @vec3=@{$_[3]};
  my @vec4;
  die "linop_parameter: 1st ($#vec1) and 3rd vector ($#vec3) have different length\n" unless ($#vec1 = $#vec3);
  for (my $i=0;$i<=$#vec1;$i++) {
    if (is_num($vec1[$i])) {
      $vec4[$i]=$vec1[$i]+$scale*($vec2[$i]-$vec3[$i]);
    } else {
      $vec4[$i]=$vec1[$i];
    }
  }
  $vec4[-1]="pending";
  $vec4[-2]=0;
  return @vec4;
}

#important
1;
