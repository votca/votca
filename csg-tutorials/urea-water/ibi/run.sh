#! /bin/bash -e

echo "Start IBI"

echo 'running csg_inverse --options "settings.xml"'
csg_inverse --options settings_IBI.xml
rm done

echo "Finished run"
# end of run file

