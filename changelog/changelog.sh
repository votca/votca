#!/bin/bash -le

shopt -s extglob

die () {
  [[ -n $1 ]] && echo "$*" >&2
  exit 1
}

[[ ${INPUT_MESSAGE} ]] || die "No message given"
[[ ${INPUT_PR_NUMBER} ]] || die "No PR number given"
message="* ${INPUT_MESSAGE#*: } (#$INPUT_PR_NUMBER)"

[[ -f CMakeLists.txt ]] || die "No CMakeLists.txt found"
version="$(sed -n 's/set(PROJECT_VERSION *"\([^"]*\)").*/\1/p' CMakeLists.txt)"
[[ ${version} ]] || die "No version found"

[[ -f CHANGELOG.md ]] || die "No CHANGELOG.md found"
version_section="$(awk "/^## Version ${version}/,/^$/{print \$0}" CHANGELOG.md)"
[[ $version_section ]] || die "Could not find section to $version"
last_line="$(echo "$version_section" | sed -n '$p')"
[[ $last_line ]] || die "Could not grep last line"

echo "Adding message '$message' after line '${last_line}'"
sed -i "/$last_line/a ${message}" CHANGELOG.md

git add CHANGELOG.md
git config --global user.name "Votca Bot"
git config --global user.email "github@votca.org"
git commit -m "Update CHANGELOG.md"
git push
