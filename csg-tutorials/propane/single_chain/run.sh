#! /usr/bin/env -S bash -e

gmx=$(type -p gmx || type -p gmx_d || { echo "gmx not found" >&2; exit 1; })

$gmx grompp -v

$gmx mdrun -v

