#! /bin/bash
#
# Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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
This script calcs the pressure for espresso and writes it to outfile

Usage: ${0##*/} outfile

Used external packages: espresso
EOF
   exit 0
fi

[ -z "$1" ] && die "${0##*/}: Missing argument"

# Espresso config file (required for certain parameters, e.g. box size)
esp="$(csg_get_property cg.inverse.espresso.blockfile "conf.esp.gz")"
[ -f "$esp" ] || die "${0##*/}: espresso blockfile '$esp' not found"

p_file="$(critical mktemp esp.pressure.val.XXXXX)"
esp_bin="$(csg_get_property cg.inverse.espresso.pressure_command "Espresso_bin")"
[ -n "$(type -p $esp_bin)" ] || die "${0##*/}: esp_bin binary '$esp_bin' not found"

esp_script="$(critical mktemp esp.pressure.tcl.XXXXX)"
esp_success="$(critical mktemp esp.pressure.done.XXXXX)"

echo "Calculating pressure"

cat > $esp_script <<EOF
set esp_in [open "|gzip -cd $esp" r]
while { [blockfile \$esp_in read auto] != "eof" } { }
close \$esp_in

set p [analyze pressure total]
set p_out [open $p_file w]
puts \$p_out $\p
close \$p_out

set out [open $esp_success w]
close \$out
EOF

critical $esp_bin $esp_script
[ -f "$esp_success" ] || die "${0##*/}: Espresso calc pressure did not end successfully. Check log."

p_now="$(cat $p_file)"

[ -z "$p_now" ] && die "${0##*/}: Could not get pressure from simulation"
echo "${p_now}" > "$1"
