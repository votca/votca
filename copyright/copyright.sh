#!/bin/bash -le

shopt -s extglob

die () {
  [[ -n $1 ]] && echo "$*" >&2
  exit 1
}

[[ ${GITHUB_BASE_REF} ]] || die "No base branch"

filelist="$(git diff --name-status --diff-filter=AMC "${GITHUB_BASE_REF}" | awk '{print $2}')"
new_date="$(date +%Y)"
echo "Updating Copyright date to ${new_date} in ${filelist}"
sed -i 's/Copyright \(.*\)-.* The VOTCA/Copyright \1-'"$new_date"' The VOTCA/' "${filelist}"

git add -u
git config user.name "Votca Bot"
git config user.email "github@votca.org"
git commit -m "Update copyright"
git push
