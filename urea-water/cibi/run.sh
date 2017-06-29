#! /bin/bash -e

echo "Running C-IBI"

echo 'running csg_inverse --options "settings.xml"'
csg_inverse --options settings.xml
