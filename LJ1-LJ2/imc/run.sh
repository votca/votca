echo 'running csg_inverse --options "settings_pre.xml"'
csg_inverse --options settings_pre.xml
rm done
echo 'running csg_inverse --options "settings.xml"'
csg_inverse --options settings.xml
