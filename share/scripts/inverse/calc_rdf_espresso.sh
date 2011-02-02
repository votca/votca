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
# Unless required by applicale law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

if [ "$1" = "--help" ]; then
    cat <<EOF
${0##*/}, version %version%
This script calcs the rdf for espresso

Usage: ${0##*/}

Used external packages: espresso
EOF
    exit 0
fi

# Topology+Trajectory read by Espresso
top_traj="$(csg_get_property cg.inverse.espresso.traj "top_traj.esp")"

# Number of snapshots before statistics are taken into account
equi_snapshots="$(csg_get_property cg.inverse.espresso.equi_snapshots 0)"

# Espresso config file (required for certain parameters, e.g. box size)
esp="$(csg_get_property cg.inverse.espresso.blockfile "conf.esp.gz")"
[ -f "$esp" ] || die "${0##*/}: espresso blockfile '$esp' not found"

esp_bin="$(csg_get_property cg.inverse.espresso.bin "Espresso_bin")" 
[ -n "$(type -p $esp_bin)" ] || die "${0##*/}: esp_bin binary '$esp_bin' not found"

type1=$(csg_get_interaction_property type1)
type2=$(csg_get_interaction_property type2)
name=$(csg_get_interaction_property name)
binsize=$(csg_get_interaction_property step)
min=$(csg_get_interaction_property min)
max=$(csg_get_interaction_property max)
index1=$(csg_get_interaction_property inverse.espresso.index1)
index2=$(csg_get_interaction_property inverse.espresso.index2)

echo "Analyzing rdf for ${type1}-${type2}"
if is_done "rdf-$name"; then
    echo "rdf analsysis for ${type1}-${type2} is already done"
else
    # Output ${name}.dist.new.tab. Calculated by Espresso.
    esp_script="$(mktemp esp.rdf.tcl.XXXXX)" 
    esp_success="$(mktemp esp.rdf.done.XXXXX)"
    cat > $esp_script <<EOF
puts "Calculating RDF. Please wait..."
# First read the original conf.esp file to get the box size
set esp_in [open "|gzip -cd $esp" r]
while { [blockfile \$esp_in read auto] != "eof" } { }
close \$esp_in


set in [open $top_traj r]
set numbins [expr int(($max-$min)/($binsize*1.))]
set bf_count 0
set rdf ""
while { [blockfile \$in read auto] != "eof" } {
  if { \$bf_count > [expr 2 + $equi_snapshots] } {
    analyze append
    # <rdf> if there's only one molecule, <rdf-intermol> otherwise
    if { \$num_molecules == 1 } {
      set rdf [analyze <rdf> [set $index1] [set $index2] $min $max \$numbins]
    } else {
      set rdf [analyze <rdf-intermol> [set $index1] [set $index2] $min $max \$numbins]
    }
  }
  incr bf_count
}
close \$in

foreach value [lindex \$rdf 1] {
  lappend rlist [lindex \$value 0]
  lappend rdflist [lindex \$value 1]
}

set out [open $name.dist.new.tab w]
foreach r \$rlist rdf \$rdflist { puts \$out "\$r \$rdf" }
close \$out

puts "Calculation finished."

set out [open $esp_success w]
close \$out
EOF
    
    critical $esp_bin $esp_script
    [ -f "$esp_success" ] || die "${0##*/}: Espresso calc rdf did not end successfully. Check log."
    
    comment="$(get_table_comment)"
    critical csg_resample --in ${name}.dist.new.tab --out ${name}.dist.new --grid ${min}:${binsize}:${max} --comment "$comment"
    mark_done "rdf-$name"
fi
