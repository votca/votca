#! /bin/bash
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

if [ "$1" = "--help" ]; then
cat <<EOF
${0##*/}, version %version%
This script runs espresso
for the Inverse Boltzmann Method

Usage: ${0##*/}

USES: run_or_exit Espresso_bin use_mpi csg_get_property check_deps

NEEDS: cg.inverse.espresso.n_steps cg.inverse.espresso.n_snapshots

OPTIONAL: cg.inverse.espresso.blockfile
EOF
   exit 0
fi

esp="$(csg_get_property cg.inverse.espresso.blockfile "conf.esp")"
[ -f "$esp" ] || die "${0##*/}: espresso blockfile '$esp' not found"

n_steps="$(csg_get_property cg.inverse.espresso.n_steps)"
[ -z "$n_steps" ] && die "${0##*/}: Could not read espresso property n_steps"

n_snapshots="$(csg_get_property cg.inverse.espresso.n_snapshots)"
[ -z "$n_snapshots" ] && die "${0##*/}: Could not read espresso property n_snapshots"


check_deps "$0"

# Topology+trajectory file
traj_esp="top_traj.esp"

# Make sure all particle indexes have been loaded into the blockfile
index_vars=$(for_all non-bonded \
csg_get_interaction_property inverse.espresso.index1)
index_vars="$index_vars $(for_all non-bonded \
csg_get_interaction_property inverse.espresso.index2)"
index_vars=$(for i in $index_vars; do echo $i; done | sort -u)
for i in $index_vars; do
    [ `grep -c $i $esp` = "1" ] || die "${0##*/}: can't find index list: $i"
done

# Calculate TCL variables: num_molecules and num_atoms and add them to
# confout.esp in case they are not already present


# load blockfile into Espresso, then integrate for $n_steps steps, then save blockfile
cat > temp_script_run_esp.tcl <<EOF
set in [open $esp r]
while { [blockfile \$in read auto] != "eof" } {}
close \$in

# Are num_molecules and num_atoms both already defined?
if { ![info exists num_molecules] || ![info exists num_atoms] } {
  # Loop over all particles to calculate num_molecules and num_atoms
  set moltypes ""
  set num_molecules 0
  set num_atoms ""
  set num_atoms_mol 0
  for { set j 0 } { \$j < [setmd n_part] } { incr j } {
    # look for molecules
    if {[lsearch \$moltypes [part \$j print molecule]] == -1} {
      lappend moltypes [part \$j print molecule]
      incr num_molecules
      if {\$j != 0 } { lappend num_atoms \$num_atoms_mol }
      set num_atoms_mol 1
    } else {
      incr num_atoms_mol
    }
  }
  lappend num_atoms \$num_atoms_mol
}

 
# Main integration loop
set pos_out [open $traj_esp w]
blockfile \$pos_out write variable box_l
blockfile \$pos_out write tclvariable num_molecules
blockfile \$pos_out write tclvariable num_atoms
close \$pos_out
for { set j 0 } { \$j < $n_snapshots } { incr j } {
  integrate $n_steps
  puts "step \$j of $n_snapshots"
  puts "  [analyze energy]"
  set pos_out [open $traj_esp a]
  if { [has_feature "MASS"] } {
    blockfile \$pos_out write particles {id type molecule mass pos v}
  } else {
    blockfile \$pos_out write particles {id type molecule pos v}
  }
  close \$pos_out
}

set out [open "confout.esp" w]
blockfile \$out write variable all
blockfile \$out write interactions
blockfile \$out write thermostat
blockfile \$out write tclvariable num_molecules
blockfile \$out write tclvariable num_atoms
blockfile \$out write tclvariable {$index_vars}
if { [has_feature "MASS"] } {
  blockfile \$out write particles {id type molecule mass pos v}
} else {
  blockfile \$out write particles {id type molecule pos v}
}
close \$out
EOF

run_or_exit Espresso_bin temp_script_run_esp.tcl
run_or_exit rm -f temp_script_run_esp.tcl
