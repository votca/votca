#!/bin/bash -le

shopt -s extglob

die () {
  [[ -n $1 ]] && echo "$*" >&2
  exit 1
}

[[ ${INPUT_MESSAGE} ]] || die "No message given"
[[ ${INPUT_PR_NUMBER} ]] || die "No PR number given"

[[ -f CMakeLists.txt ]] || die "No CMakeLists.txt found"
version="$(sed -n 's/set(PROJECT_VERSION *"\([^"]*\)").*/\1/p' CMakeLists.txt)"
[[ ${version} ]] || die "No version found"

CHANGELOG=CHANGELOG.rst
version_section="$(awk -v r="^Version ${version}( |$)" '($0 ~ "^Version"){go=0} ($0 ~ r){go=1}{if(go==1){print $0}}' "${CHANGELOG}")"
line_nr="$(sed -n "/^Version ${version}\( \|$\)/=" "${CHANGELOG}")"
[[ -z ${version_section} ]] || break
message="-  ${INPUT_MESSAGE#*: } (#$INPUT_PR_NUMBER)"

[[ $version_section ]] || die "Could not find section to $version"
echo "Found section for $version (starting line ${line_nr})"
last_line="$(echo "$version_section" | sed '/^[[:space:]]*$/d' | sed -n '$p')"
[[ $last_line ]] || die "Could not grep last line"

if [[ ${last_line} = -* || ${last_line} = '   '[^\ ]* ]]; then
  echo "Adding message '$message' after line '${last_line}'"
  sed -i "/$last_line/a ${message}" "${CHANGELOG}"
elif [[ -z ${last_line//=} ]]; then #section header
  (( line_nr=line_nr+2 ))
  echo "Adding message '$message' after line ${line_nr}"
  sed -i "${line_nr}a ${message}\n" "${CHANGELOG}"
else
  die "Last line isn't an item (does not start with -) nor a section header (===), case not implemented"
fi

[[ ${CI} = 'true' ]] || exit 0
git add "${CHANGELOG}"
git config user.name "Votca Bot"
git config user.email "github@votca.org"
git commit -m "Update ${CHANGELOG}"
git push
