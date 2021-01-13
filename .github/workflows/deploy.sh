#!/bin/bash -e

builddir="$PWD"

pushd "devdoc"
rm -rf -- *
git checkout -- CNAME manual.pdf
cp -r "${builddir}/share/doc/html" .
mv html/* .
rmdir html
git add --all .
git config user.name "Votca Bot"
git config user.email "github@votca.org"
git commit -m "Documentation Update for votca/votca@${GITHUB_SHA}" || true
if [[ ${GITHUB_REF} = "refs/heads/master" ]]; then
  git push origin gh-pages
else
  git show --no-color | head -c 1k
fi
popd

if [[ ${GITHUB_BASE_REF} = stable || ${GITHUB_REF} = refs/heads/stable || ${BRANCH} = stable ]]; then
  exit 0
fi

pushd "userdoc"
rm -rf -- * .buildinfo .doctrees
git checkout -- CNAME
cp -r "${builddir}/sphinx.html" .
mv sphinx.html/* sphinx.html/.buildinfo sphinx.html/.doctrees sphinx.html/.nojekyll .
rmdir sphinx.html
git add --all .
git config user.name "Votca Bot"
git config user.email "github@votca.org"
git commit -m "Documentation Update for votca/votca@${GITHUB_SHA}" || true
if [[ ${GITHUB_REF} = "refs/heads/master" ]]; then
  git push origin master
else
  git show --no-color | head -c 1k
fi
popd
