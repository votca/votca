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
  try_versions=( ${version%-dev}-rc.1 )
elif [[ $version = 20??.* ]]; then
  try_versions=( "${version%.*}.$((${version##*.}+1))" )
elif [[ $version = 20?? ]]; then
  try_versions=( "${version}.1" )
elif [[ $version = *-rc.* ]]; then
  try_versions=( "${version%-rc*}-rc.$((${version#*-rc.}+1))" )
else
  die "Unknown version scheme, found $version"
fi
try_versions+=( $version )
echo "Trying versions ${try_versions[@]}"

if [[ -f CHANGELOG.rst ]]; then
  CHANGELOG=CHANGELOG.rst
  for v in ${try_versions[@]}; do
    version_section="$(awk -v r="^Version ${v}( |$)" '($0 ~ "^Version"){go=0} ($0 ~ r){go=1}{if(go==1){print $0}}' "${CHANGELOG}")"
    [[ -z ${version_section} ]] || break
  done
  message="-  ${INPUT_MESSAGE#*: } (#$INPUT_PR_NUMBER)"
else
  die "No supported CHANGELOG found"
fi

[[ $version_section ]] || die "Could not find section to $version"
echo "Found section for $v"
last_line="$(echo "$version_section" | sed '/^[[:space:]]*$/d' | sed -n '$p')"
[[ $last_line ]] || die "Could not grep last line"
[[ -z ${last_line##-*} ]] || die "Last line isn't an item (does not start with -)"

echo "Adding message '$message' after line '${last_line}'"
sed -i "/$last_line/a ${message}" "${CHANGELOG}"

[[ ${CI} = 'true' ]] || exit 0
git add "${CHANGELOG}"
git config user.name "Votca Bot"
git config user.email "github@votca.org"
git commit -m "Update ${CHANGELOG}"
git push
