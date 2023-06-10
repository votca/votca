#!/bin/bash -e
set -euo pipefail

echo 'running the "pre" iteration, generating dcdh.npz with a longer cut-off in step_000'
pushd pre
csg_inverse --options settings.xml
popd

echo 'copying files from pre'
for f in pre/step_000/*.pot.new; do
    fname="$(basename "$f")"
    dest="${fname%%new}in"
    cp $f $dest
    echo "copied $f to $dest"
done
cp pre/dcdh.npz .
echo "copied pre/dcdh.npz to dcdh.npz"

echo 'running the main iterations with shorter cut-off (faster)'
csg_inverse --options settings.xml
