md5sum *.pl *.sh linsolve.m linsolve.octave csg_table > MD5SUM
sed -ie '/[:space:]*recalc_md5.sh$/d' MD5SUM

