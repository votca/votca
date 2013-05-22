#!/bin/bash -e

VOTCASHARE="$(csg_call --show-share)"

texfile="$PWD/calculators.tex"
rm -f $texfile; touch $texfile

for package in ctp kmc; do

	calculators="$(${package}_run --list | sed -ne 's/^\s\+\([a-z]*\)\s*\(.*\)/\1/p')"

	# loop over all calculators
	for calculator in ${calculators}; do

		xmlfile=${VOTCASHARE}/${package}/xml/$calculator.xml

		if [ ! -f "$xmlfile" ]; then 
			continue
		fi

		votca_property --file $xmlfile --format TEX >> $texfile

	done
done
