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
This script implemtents the potential update for the PMF method

Usage: ${0##*/}

USES: do_external for_all csg_get_interaction_property check_deps csg_get_property

NEEDS: type inter cg.inverse.espresso.inclusion_type

OPTIONAL: cg.inverse.espresso.analysis cg.inverse.program
EOF
   exit 0
fi

check_deps "$0"

# pmf analysis input
err_ptypes="$(csg_get_property cg.inverse.espresso.analysis "err.ptypes.dat")"
[ -f "$err_ptypes" ] || die "${0##*/}: PMF analysis file $err_ptypes not found"
err_ptypes="$(cat $err_ptypes)"

sim_prog="$(csg_get_property cg.inverse.program)"
# convert to lower case
sim_prog="$(echo $sim_prog | tr [:upper:] [:lower:])"

# particle types involved in the optimization (except inclusion)
ptypes="$(for_all non-bonded csg_get_interaction_property type)"

update_factor="$(csg_get_property cg.inverse.espresso.update_factor 0.5)"
[ -z "$update_factor" ] && die "${0##*/}: Could not read metadynamics property update_factor"


# Inclusion
incl_type="$(csg_get_property cg.inverse.espresso.inclusion_type)"

# Check that err_ptypes and ptypes have the same number of elements
if [ "$(echo $err_ptypes | wc -w)" -ne "$(echo $ptypes | wc -w)" ]; then
    die "Error in the PMF analysis file"
fi

if [ "$sim_prog" = "espresso" ]; then
    
    # Espresso config file (required for certain parameters, e.g. box size)
    esp="confout.esp.gz"
    [ -f "$esp" ] || die "${0##*/}: espresso blockfile '$esp' not found"

    esp_bin="$(csg_get_property cg.inverse.espresso.bin "Espresso_bin")"

    esp_script="$(mktemp esp.pmf_pot.tcl.XXXXX)"
    esp_success="$(mktemp esp.pmf_pot.done.XXXXX)"
    cat > $esp_script <<EOF
puts "Updating potentials."
# First read the original conf.esp file to get the box size
set esp_in [open "|gzip -cd $esp" r]
while { [blockfile \$esp_in read auto] != "eof" } { }
close \$esp_in


for { set i 0 } { \$i < [llength {$ptypes}] } { incr i } {
  set interaction [inter [lindex {$ptypes} \$i] $incl_type]
  # Update epsilon of LJ-type interaction
  set parameter [lindex \$interaction 3]
  # Apply update factor
  set parameter [expr \$parameter * ($update_factor * [lindex {$err_ptypes} \$i] + 1-$update_factor)]
  lset interaction 3 \$parameter
  eval [ concat inter \$interaction ]

}

# Update blockfile
set out [open "|gzip -c - > confout.esp.gz" w]
blockfile \$out write variable all
blockfile \$out write interactions
blockfile \$out write thermostat
blockfile \$out write tclvariable num_molecules
blockfile \$out write tclvariable num_atoms
if { [has_feature "MASS"] } {
  if { [has_feature "VIRTUAL_SITES"] } {
    blockfile \$out write particles {id type molecule mass virtual pos v}
  } else {
    blockfile \$out write particles {id type molecule mass pos v}
  }
} else {
  if { [has_feature "VIRTUAL_SITES"] } {
    blockfile \$out write particles {id type molecule virtual pos v}
  } else {
    blockfile \$out write particles {id type molecule pos v}
  }
}
close \$out


set out [open $esp_success w]
puts \$out "done"
close \$out
EOF

    run_or_exit $esp_bin $esp_script
    [ -f "$esp_success" ] || die "${0##*/}: Espresso calc rdf did not end successfully. Check log."    

else
    die "PMF analysis only supported on ESPResSo"
fi
