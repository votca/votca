#! /bin/bash 

not_to_check="inverse.sh inverse_start.sh errorbars.sh ${0##*/} functions_gromacs.sh functions_inverse.sh"
[[ -z "$1" ]] && echo help with ${0##*/} --help && exit 
if [ "$1" = "--help" ]; then
  echo Usage: ${0##*/} WORD FILE1 FILE2 ...
  echo Checks if a program is in USES block of a file
  echo "Example ${0##*/} awk *.sh"
  echo WORD = -- does some magic
  echo "WORD is set to \$(for i in *.sh; do ./\$i --help 2>&1 | sed -n 's/USES: \+\(.*\) *$/\1/p' | sed 's/ /\n/g'; done | sort | uniq)"
  echo It will always ignore: $not_to_check
  exit 0
fi

if [ "$1" = "--" ]; then
  whates="$(for i in *.sh; do [[ -z "${not_to_check##*$i*}" ]] && continue; ./$i --help 2>&1 | sed -n 's/USES: \+\(.*\) *$/\1/p' | sed 's/ /\n/g'; done | sort | uniq )"
else
  whates="$1"
fi
shift

echo what: $whates
for i in $@; do
  [[ -z "${not_to_check##*$i*}" ]] && continue
  echo Checking $i
  [[ ! -x "$i" ]] && echo "$i is not executable" && continue
  ./$1 --help &> /dev/null || { echo "$i has no help"; continue; }
  [[ -z "(./$1 --help | grep "USES:")" ]] && echo "$i has no USES in help" && continue
  for what in $whates; do
    [[ -z "${what##\$*}" ]] && what="\\$what"
    #what found in file and uses -> ok
    [[ -n "$(grep -v "USES:" "$i" | grep -Ee "$what([[:space:]]|\"|$)")" ]] && \
      [[ -n "$(./$i --help | grep -Ee "USES:.*$what([[:space:]]|$)")" ]] && \
      continue
    #what found in file, but not in uses
    [[ -n "$(grep -v "USES:" "$i" | grep -Ee "$what([[:space:]]|\"|$)")" ]] && \
      [[ -z "$(./$i --help | grep -Ee "USES:.*$what([[:space:]]|$)")" ]] && \
      echo "$i: $what found, but NOT in USES -> add it"
    #what not found in file, but in uses
    [[ -z "$(grep -v "USES:" "$i" | grep -Ee "$what([[:space:]]|\"|$)")" ]] && \
      [[ -n "$(./$i --help | grep -Ee "USES:.*$what([[:space:]]|$)")" ]] && \
      echo "$i: $what found in USES, but NOT in content -> remove it"
    #what not found in file and uses -> ok
    [[ -z "$(grep -v "USES:" "$i" | grep -Ee "$what([[:space:]]|\"|$)")" ]] && \
      [[ -z "$(./$i --help | grep -Ee "USES:.*$what([[:space:]]|$)")" ]] && \
      continue
  done
done
