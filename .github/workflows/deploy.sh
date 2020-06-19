#!/bin/bash -e

builddir="$PWD/builddir"

git clone --depth=1 https://github.com/votca/doxygen.git "$HOME/devdoc"
pushd "$HOME/devdoc"
rm -rf -- *
git checkout -- CNAME keys
cp -r "${builddir}/share/doc/html" .
mv html/* .
rmdir html
cp "${builddir}/csg-manual/manual.pdf" .
git add --all .
git config --global user.name "Votca Bot"
git config --global user.email "github@votca.org"
git commit -m "Documentation Update for votca/votca@${GITHUB_SHA}"
if [[ ${GITHUB_REF} = "refs/heads/master" && ${VOTCA_BOT_TOKEN} ]]; then
  git push "https://${VOTCA_BOT_TOKEN}@github.com/votca/doxygen.git" gh-pages
else
  git show --no-color | head -c 1k
fi
popd

git clone --depth=1 https://github.com/votca/votca.github.io.git "$HOME/userdoc"
pushd "$HOME/userdoc"
rm -rf -- * .buildinfo .doctrees
git checkout -- deploy.enc
cp -r "${builddir}/sphinx.html" .
mv sphinx.html/* sphinx.html/.buildinfo sphinx.html/.doctrees sphinx.html/.nojekyll .
rmdir sphinx.html
git add --all .
git config --global user.name "Votca Bot"
git config --global user.email "github@votca.org"
git commit -m "Documentation Update for votca/votca@${GITHUB_SHA}"
if [[ ${GITHUB_REF} = "refs/heads/master" && ${VOTCA_BOT_TOKEN} ]]; then
  git push "https://${VOTCA_BOT_TOKEN}@github.com/votca/votca.github.io.git" master
else
  git show --no-color | head -c 1k
fi
popd
