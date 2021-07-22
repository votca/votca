#!/bin/bash

set -e

# pre-iterations with IBI
cmd='csg_inverse --options "settings-ibi.xml"'
echo "now running: ${cmd}"
$cmd
rm done

# IMC iterations
cmd='csg_inverse --options "settings-imc.xml"'
echo "now running: ${cmd}"
$cmd
