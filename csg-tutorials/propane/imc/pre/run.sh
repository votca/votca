#!/usr/bin/env bash
set -e

# pre-iterations with IBI
cmd='csg_inverse --options settings.xml'
echo "now running: ${cmd}"
$cmd

echo "now copying final potentials from IBI for IMC"
for nb in A-A B-B A-B; do
    if [[ ! -f ../${nb}.pot.in ]]; then
        cp "step_005/${nb}.pot.new" "../${nb}.pot.in"
    fi
done
