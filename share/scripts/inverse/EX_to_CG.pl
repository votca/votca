#! /usr/bin/perl -w

use strict;

die "2 more parameters are nessary\n" if ($#ARGV<1);

my $file="$ARGV[0]";
open(FILE1,$file) or die "$file not found\n";

my $file2="$ARGV[1]";
open(FILE2,"> $file2") or die "$file2 not found\n";
<FILE1>;
print FILE2 "CG Water\n";
my $n_atoms=<FILE1>;
chomp($n_atoms);
print "$n_atoms Atoms found\n";
my $n_mols=$n_atoms/3;
print FILE2 " $n_mols\n";

my %mass;
$mass{OW}=15.9994;
$mass{HW1}=1.008;
$mass{HW2}=1.008;
$mass{Cl}=35.45300;
$mass{Na}=22.98977;

my %cg_name;
$cg_name{OW}="CG";
$cg_name{HW1}="CG";
$cg_name{HW2}="CG";
$cg_name{Cl}="Cl";
$cg_name{Na}="Na";

my $this_mol_nr;
my $this_mol_name;
my $this_atom_name;

my $mol_nr=0;
my $mol_name;
my $atom_name;

my @cg_vel;
my @cg_pos;
my $cg_mass;

for (my $atom=1;$atom<=$n_atoms;$atom++) {
   my @lineparts=split(' ',<FILE1>);
   #print $#lineparts."\n";
   die "Atom is not supported\n" unless ($lineparts[0] =~ /^([0-9]*)(SOL|Cl|Na)$/);
   $this_mol_nr=$1;
   $this_mol_name=$2;
   $this_atom_name=$lineparts[1];
   die "Unknown Atom type $this_atom_name\n" unless ($mass{$this_atom_name});
   #new mole write out !!!

   if ($this_mol_nr != $mol_nr ){
      if ( $mol_nr != 0) {
         foreach (0..2) {
            $cg_pos[$_]/=$cg_mass;
            $cg_vel[$_]/=$cg_mass;
         }
         printf FILE2 "%5i%-5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
               $mol_nr,$mol_name,$cg_name{$atom_name},$mol_nr,$cg_pos[0],$cg_pos[1],$cg_pos[2],$cg_vel[0],$cg_vel[1],$cg_vel[2];
      }
      $mol_nr=$this_mol_nr;
      $mol_name=$this_mol_name;
      $atom_name=$this_atom_name;
      $cg_mass=0.0;
      foreach (0..2) {
         $cg_pos[$_]=0.0;
         $cg_vel[$_]=0.0;
      }
   }
   $cg_mass+=$mass{$lineparts[1]};
   foreach (0..2) {
      $cg_pos[$_]+=$mass{$lineparts[1]}*$lineparts[$_+3];
      $cg_vel[$_]+=$mass{$lineparts[1]}*$lineparts[$_+6];
   }
}
#print last mol
foreach (0..2) {
   $cg_pos[$_]/=$cg_mass;
   $cg_vel[$_]/=$cg_mass;
}
printf FILE2 "%5i%-5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
              $mol_nr,$mol_name,$cg_name{$atom_name},$mol_nr,$cg_pos[0],$cg_pos[1],$cg_pos[2],$cg_vel[0],$cg_vel[1],$cg_vel[2];

print FILE2 <FILE1>;
close(FILE1) or die "Error at closing $file\n";
close(FILE2) or die "Error at closing $file2\n";