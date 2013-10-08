#! /bin/bash -e

echo "Start KB-IBI"

csg_inverse --options settings_KB-IBI_1.xml
rm done
csg_inverse --options settings_KB-IBI_2.xml
rm done
csg_inverse --options settings_KB-IBI_3.xml
rm done

echo "Finished run"
# end of run file

