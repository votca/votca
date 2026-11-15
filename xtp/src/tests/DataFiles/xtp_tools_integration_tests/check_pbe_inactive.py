#!/usr/bin/env python3
#
# Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""Check GW-BSE results: PbE with inactive screening, streamed RPA."""

import argparse
import sys

import h5py
import numpy as np


def read(path):
    """Return QP energies by level index, the HOMO and the singlets."""
    with h5py.File(path, "r") as f:
        q = f["QMdata"]
        qpmin = int(q.attrs["qpmin"][0])
        homo = int(q.attrs["number_alpha_electrons"][0]) - 1
        qp = q["QPpert_energies"][()].ravel()
        singlets = q["BSE_singlet"]["eigenvalues"][()].ravel()
    return {qpmin + i: e for i, e in enumerate(qp)}, homo, singlets


def main():
    """Compare a PbE run with the full one, or check that it is sane."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pbe", help="orb file of the PbE run")
    parser.add_argument("--full", help="orb file of the full run: with "
                        "exact screening and complete G the virtual QP "
                        "energies must equal the full ones")
    parser.add_argument("--same", help="orb file of a run that must give "
                        "the same QP and singlet energies (all levels)")
    parser.add_argument("--atol", type=float, default=2e-4,
                        help="tolerance in Hartree")
    args = parser.parse_args()

    qp, homo, singlets = read(args.pbe)
    values = np.array(list(qp.values()))
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(singlets)):
        print("non-finite QP or BSE energies")
        return 1
    if len(singlets) == 0 or singlets[0] <= 0.0:
        print("no positive singlet")
        return 1
    if qp[homo + 1] <= qp[homo]:
        print("LUMO below HOMO")
        return 1
    if args.same:
        other, other_homo, other_singlets = read(args.same)
        if other_homo != homo or set(other) != set(qp):
            print("different levels")
            return 1
        diff = max(abs(qp[k] - other[k]) for k in qp)
        n = min(len(singlets), len(other_singlets))
        sdiff = float(np.max(np.abs(singlets[:n] - other_singlets[:n])))
        print(f"{len(qp)} QP levels, max deviation {diff:.2e} Ha; "
              f"{n} singlets, max deviation {sdiff:.2e} Ha")
        if diff > args.atol or sdiff > args.atol:
            return 1
    if args.full:
        full, full_homo, _ = read(args.full)
        if full_homo != homo:
            print(f"HOMO index differs: {homo} vs {full_homo}")
            return 1
        virt = [k for k in qp if k > homo and k in full]
        diff = max(abs(qp[k] - full[k]) for k in virt)
        print(f"{len(virt)} virtual levels, max deviation {diff:.2e} Ha")
        if diff > args.atol:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
