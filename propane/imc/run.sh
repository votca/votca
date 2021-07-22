#!/bin/bash

set -e

ln -sf grompp-ibi.mdp grompp.mdp
echo 'running csg_inverse --options "settings_pre.xml"'
csg_inverse --options settings-preibi.xml
rm done
ln -sf grompp-imc.mdp grompp.mdp
echo 'running csg_inverse --options "settings.xml"'
csg_inverse --options settings.xml
rm -f grompp.mdp
