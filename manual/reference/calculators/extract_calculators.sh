#!/bin/bash -e

texfile="$PWD/calculators.tex"
rm -f $texfile; touch $texfile

for package in xtp_tools xtp_run xtp_parallel xtp_dump; do
	exe=$(echo "$@" | xargs -n1 echo | grep "$package")
	calculators="$(${exe} --list | sed -ne 's/^\s\+\([a-z,0-9]*\)\s*\(.*\)/\1/p')"
	echo $calculators

	# loop over all calculators
	for calculator in ${calculators}; do
		xmlfile=../../../share/xml/${calculator}.xml
		echo $calculator
		echo $xmlfile
		if [ ! -f "$xmlfile" ]; then 
			continue
		fi
		echo "votca_property --file $xmlfile --format TEX --level 2"
		votca_property --file $xmlfile --format TEX --level 2 >> $texfile

	done
done
