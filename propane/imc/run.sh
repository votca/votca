#!/bin/bash

set -e

ln -sf settings-ibi.mdp settings.xml
ln -sf grompp-ibi.mdp grompp.mdp
echo 'running csg_inverse --options "settings.xml"'
csg_inverse --options settings.xml
rm done

ln -sf settings-imc.mdp settings.xml
ln -sf grompp-imc.mdp grompp.mdp
echo 'running csg_inverse --options "settings.xml"'
csg_inverse --options settings.xml
rm -f settings.xml grompp.mdp
