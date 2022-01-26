#!/bin/bash -e
set -euo pipefail

echo 'running the "pre1-step0" iteration, generating potential guess'
pushd pre1-step0
! csg_inverse --options settings.xml  # it will fail after step_000 due to iterations_max=-1; exclamation mark prevents stop
popd

echo 'copying files from pre1-step0/step_000 to pre2-ibi'
for f in pre1-step0/step_000/*.pot.new; do
    fname="$(basename "$f")"
    dest="pre2-ibi/${fname%%new}in"
    cp $f $dest
    echo "copied $f to $dest"
done

echo 'running the "pre2-ibi" iterations. Running hncn directly would lead to instability.'
pushd pre2-ibi
csg_inverse --options settings.xml
popd

echo 'copying files from pre2-ibi/step_004 to main dir'
for f in pre2-ibi/step_004/*.pot.new; do
    fname="$(basename "$f")"
    dest="${fname%%new}in"
    cp $f $dest
    echo "copied $f to $dest"
done

echo 'running the main iterations'
csg_inverse --options settings.xml
