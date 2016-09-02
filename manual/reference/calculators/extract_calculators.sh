a#!/bin/bash -e

VOTCASHARE="$(csg_call --show-share)"

texfile="$PWD/calculators.tex"
rm -f $texfile; touch $texfile

for package in xtp_tools xtp_run xtp_parallel xtp_dump xtp_kmc_run; do
	calculators="$(${package} --list | sed -ne 's/^\s\+\([a-z,0-9]*\)\s*\(.*\)/\1/p')"
	echo $calculators

	# loop over all calculators
	for calculator in ${calculators}; do
		library="$(echo ${package} | sed -e 's/\_[^\_]*$//g')"
		xmlfile=${VOTCASHARE}/${library}/xml/$calculator.xml
		echo $calculator
		echo $xmlfile
		if [ ! -f "$xmlfile" ]; then 
			continue
		fi
		echo "votca_property --file $xmlfile --format TEX --level 2"
		votca_property --file $xmlfile --format TEX --level 2 >> $texfile

	done
done
