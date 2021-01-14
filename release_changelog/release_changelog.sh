#!/bin/bash -le

shopt -s extglob

die () {
  [[ -z $1 ]] || echo "$*" >&2
  exit 1
}

if [[ -z ${INPUT_VERSION} ]]; then
  [[ -f CMakeLists.txt ]] || die "No CMakeLists.txt found"
  version="$(sed -n 's/set(PROJECT_VERSION *"\([^"]*\)").*/\1/p' CMakeLists.txt)"
  [[ ${version} ]] || die "No version found"
else
  version="${INPUT_VERSION#refs/tags/v}"
fi
[[ -n ${INPUT_CHANGELOG} ]] || die "INPUT_CHANGELOG is empty"

if [[ -f CHANGELOG.rst ]]; then
  CHANGELOG=CHANGELOG.rst
  version="${version//_/\\\\\\\\_}" # backslash underscores
  version_section="$(awk -v r="^Version ${version}( |$)" '($0 ~ "^Version"){go=0} ($0 ~ r){go=1}{if(go==1){print $0}}' "${CHANGELOG}" | sed -e '1,2d' -e '/^$/d')"
else
  die "No supported CHANGELOG found"
fi

[[ $version_section ]] || die "Could not find section to $version"
echo "Writing following changelog to $INPUT_CHANGELOG:"
echo "$version_section"
echo "$version_section" > "${INPUT_CHANGELOG}"
