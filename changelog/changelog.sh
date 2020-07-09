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

if [[ -f CHANGELOG.rst ]]; then
  CHANGELOG=CHANGELOG.rst
  version="${version//_/\\\\\\\\_}" # backslash underscores
  version_section="$(awk -v r="^Version ${version}( |$)" '($0 ~ "^Version"){go=0} ($0 ~ r){go=1}{if(go==1){print $0}}' "${CHANGELOG}")"
  message="-  ${INPUT_MESSAGE#*: } (#$INPUT_PR_NUMBER)"
elif [[ -f CHANGELOG.md ]]; then
  CHANGELOG=CHANGELOG.md
  version_section="$(awk "/^## Version ${version}( |$)/,/^$/{print \$0}" "${CHANGELOG}")"
  message="* ${INPUT_MESSAGE#*: } (#$INPUT_PR_NUMBER)"
else
  die "No supported CHANGELOG found"
fi

[[ $version_section ]] || die "Could not find section to $version"
last_line="$(echo "$version_section" | sed -n '$p')"
[[ $last_line ]] || die "Could not grep last line"

echo "Adding message '$message' after line '${last_line}'"
sed -i "/$last_line/a ${message}" "${CHANGELOG}"

git add "${CHANGELOG}"
git config user.name "Votca Bot"
git config user.email "github@votca.org"
git commit -m "Update ${CHANGELOG}"
git push
