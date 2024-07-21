#!/usr/bin/env -S bash -e
set -e

# IMC iterations
cmd='csg_inverse --options settings.xml'
echo "now running: ${cmd}"
$cmd
