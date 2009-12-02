#! /bin/bash 

#note the space at the beginning and the end !!!
not_to_check=" ${0##*/} "

if [ "$1" = "--debug" ]; then
  debug="yes"
  shift
else
  debug="no"
fi

if [ -z "$1" ]; then
  echo MIssing argument >&2 
  echo help with ${0##*/} --help >&2
  exit 1
fi

if [ "$1" = "--help" ]; then
  echo Usage: ${0##*/} WORD FILE1 FILE2 ...
  echo Checks if a program is in USES block of a file
  echo "Example ${0##*/} awk *.sh"
  echo WORD = -- does some magic
  echo "WORD is set to \$(for i in *.sh *.pl; do ./\$i --help 2>&1 | sed -n 's/USES: \+\(.*\) *$/\1/p' | sed 's/ /\n/g'; done | sort | uniq)"
  echo "missing FILE means '*.sh *.pl'"
  echo "So if you are LAZY just run '${0##*/} --'"
  echo It will always ignore: $not_to_check
  exit 0
fi

if [ "$1" = "--" ]; then
  whates="$(for i in *.sh *.pl; do [[ -z "${not_to_check##* $i *}" ]] && continue; ./$i --help 2>&1 | sed -n 's/\(USES:\|PROVIDES:\) \+\(.*\) *$/\2/p' | sed 's/ /\n/g'; done | sort | uniq )"
else
  whates="$1"
fi
shift

if [ -z "$1" ]; then
  set -- *.pl *.sh
fi

echo files to check: $@
echo
echo what to check: $whates
echo
for i in $@; do
  [[ -z "${not_to_check##* $i *}" ]] && continue
  echo Checking $i
  [[ ! -x "$i" ]] && echo "$i is not executable" && continue
  ./$i --help &> /dev/null || { echo "$i has no help"; continue; }
  [[ -z "(./$i --help | grep "USES:")" ]] && echo "$i has no USES in help" && continue
  for what in $whates; do
    #buildins in perl
    if [ -z "${i%%*.pl}" ]; then
      [[ "$what" = "die" ]] && continue
      [[ "$what" = "log" ]] && continue
      [[ "$what" = "printf" ]] && continue
    fi
    #variable
    if [ -z "${what##\$*}" ]; then
      if [ -z "${i%%*.sh}" ]; then
	what="${what}"
	what2="$what"
      elif [ -z "${i%%*.pl}" ]; then 
	what2="${what}"
	what="\$ENV\{'?${what##\$}'?\}"
      fi
    else
      what2="$what"
    fi
    #pattern in the content of the file
    pattern1="(^|[^a-zA-Z])$what([^a-zA-Z._]|$)"
    #pattern in the help
    pattern2="(USES:|PROVIDES:).*[[:space:]]$what2([[:space:]]|$)"
    in_help="no"
    in_content="no"
    [[ -n "$(grep -Ev "(USES:|PROVIDES:)" "$i" | grep -Ee "$pattern1")" ]] && in_content="yes"
    [[ -n "$(./$i --help | grep -Ee "$pattern2")" ]] && in_help="yes"
    if [ "$debug" = "yes" ]; then
      echo "cont $in_content" help "$in_help"
      echo "p1 $pattern1 p2 $pattern2"
    fi
    #what found in file and uses -> ok
    [[ "$in_help" = "yes" ]] && [[ "$in_content" = "yes" ]] && continue
    #what found in file, but not in uses
    [[ "$in_help" = "no" ]] && [[ "$in_content" = "yes" ]] && \
      echo "$i: $what found, but NOT in USES -> add it"
    #what not found in file, but in uses
    [[ "$in_help" = "yes" ]] && [[ "$in_content" = "no" ]] && \
      echo "$i: $what found in USES, but NOT in content -> remove it"
    #what not found in file and uses -> ok
    [[ "$in_help" = "no" ]] && [[ "$in_content" = "no" ]] && continue
  done
done
