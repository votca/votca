#!/bin/bash -e

texfile="$PWD/calculators.tex"
rm -f $texfile; touch $texfile

die() {
  echo "$*" >&2
  exit 1
}
[[ -z ${VOTCASHARE} ]] && die "${0##*/}: VOTCASHARE not set"
[[ -z ${VOTCA_PROPERTY} ]] && die "${0##*/}: VOTCA_PROPERTY not set"
export VOTCASHARE

for package in xtp_tools xtp_run xtp_parallel xtp_dump; do
	exe=$(echo "$@" | xargs -n1 echo | grep "$package")
	calculators="$(${exe} --list | sed -ne 's/^\s\+\([a-z,0-9]*\)\s*\(.*\)/\1/p')"
	echo Found calculators: $calculators

	# loop over all calculators
	for calculator in ${calculators}; do
		xmlfile=${VOTCASHARE}/xml/${calculator}.xml
		echo Using $xmlfile for $calculator
		if [ ! -f "$xmlfile" ]; then 
			continue
		fi
		echo "${VOTCA_PROPERTY} --file $xmlfile --format TEX --level 2"
		${VOTCA_PROPERTY} --file $xmlfile --format TEX --level 2 >> $texfile

	done
done
