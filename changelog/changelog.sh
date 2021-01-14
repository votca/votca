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

if [[ $version = *-dev ]]; then
  :
elif [[ $version = *.*.* ]]; then
  try_versions=( "${version%.*}.$((${version##*.}+1))" )
elif [[ $version = *.* ]]; then
  try_versions=( "${version}.1" )
elif [[ $version = *_rc* ]]; then
  try_versions=( "${version%_rc*}_rc$((${version#*_rc}+1))" )
else
  die "Unknown version scheme, found $version"
fi
try_versions+=( $version )
echo "Trying versions ${try_versions[@]}"

if [[ -f CHANGELOG.rst ]]; then
  CHANGELOG=CHANGELOG.rst
  for v in ${try_versions[@]}; do
    vb="${v//_/\\\\\\\\_}" # backslash underscores
    version_section="$(awk -v r="^Version ${vb}( |$)" '($0 ~ "^Version"){go=0} ($0 ~ r){go=1}{if(go==1){print $0}}' "${CHANGELOG}")"
    [[ ${version_section} ]] || break
  done
  message="-  ${INPUT_MESSAGE#*: } (#$INPUT_PR_NUMBER)"
else
  die "No supported CHANGELOG found"
fi

[[ $version_section ]] || die "Could not find section to $version"
echo "Found section for $v"
last_line="$(echo "$version_section" | sed '/^[[:space:]]*$/d' | sed -n '$p')"
[[ $last_line ]] || die "Could not grep last line"

echo "Adding message '$message' after line '${last_line}'"
sed -i "/$last_line/a ${message}" "${CHANGELOG}"

git add "${CHANGELOG}"
git config user.name "Votca Bot"
git config user.email "github@votca.org"
git commit -m "Update ${CHANGELOG}"
git push
