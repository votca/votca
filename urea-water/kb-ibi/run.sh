#! /bin/bash -e

echo "Start KB-IBI"

csg_inverse --do-iterations 25 --options settings.xml
csg_inverse --do-iterations 35 --options settings_2.xml
csg_inverse --options settings_3.xml

echo "Finished run"
# end of run file

