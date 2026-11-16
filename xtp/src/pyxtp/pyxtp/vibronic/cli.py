"""Subcommands of the xtp_vibronic script (the parser is in the script,
so that its help works without pyxtp).

Subcommands:

  hessian    run (or continue) the displaced force calculations; several
             processes or array-job parts can share one workdir
  status     progress of a workdir
  modes      collect the Hessian, print the frequencies, save modes.npz
  lineshape  couplings and line shape (optionally an ensemble spectrum)
"""

from __future__ import annotations

import operator
import sys
from typing import List, Optional

import numpy as np

_OPS = {">": operator.gt, "<": operator.lt, ">=": operator.ge,
        "<=": operator.le}


def _apply_options_file(options, filename):
    """Set a dftgwbse options file (<options><dftgwbse>...) on Options."""
    import xmltodict

    with open(filename) as f:
        data = xmltodict.parse(f.read())
    try:
        tree = data["options"]["dftgwbse"]
    except (KeyError, TypeError):
        sys.exit(f"{filename}: expected <options><dftgwbse>...</dftgwbse>"
                 "</options>")

    def walk(node, obj, path):
        for key, value in node.items():
            if key.startswith("@") or key.startswith("#"):
                continue
            if not hasattr(obj, key):
                sys.exit(f"{filename}: unknown option "
                         f"{'.'.join(path + [key])}")
            if isinstance(value, dict):
                walk(value, getattr(obj, key), path + [key])
            else:
                setattr(obj, key, "" if value is None else value)

    walk(tree, options, [])


def _settings(pairs: Optional[List[str]]):
    out = {}
    for p in pairs or []:
        if "=" not in p:
            sys.exit(f"--set expects key=value, got {p}")
        k, v = p.split("=", 1)
        out[k.strip().replace(".", "/")] = v.strip()
    return out


def _parse_where(expr: str):
    for op in (">=", "<=", ">", "<"):
        if op in expr:
            col, val = expr.split(op)
            return int(col), _OPS[op], float(val)
    sys.exit(f"--where: cannot parse {expr} (e.g. '1>4.0')")


def cmd_hessian(args) -> None:
    from ase.io import read

    from ..xtp import xtp
    from .hessian import FiniteDifferenceHessian, xtp_ground_state_factory

    options = None
    if args.options:
        options = xtp().options
        _apply_options_file(options, args.options)
    factory = xtp_ground_state_factory(_settings(args.set), args.threads,
                                       warm_start=not args.no_warm_start,
                                       options=options)
    if args.xyz:
        indices = None if args.atoms is None else \
            [int(i) for i in args.atoms.split(",")]
        fd = FiniteDifferenceHessian(read(args.xyz), factory, args.workdir,
                                     args.step, args.scheme, indices)
    else:
        fd = FiniteDifferenceHessian.open(args.workdir, factory)
    n = len(fd.displacements())
    print(f"{n} force calculations in {fd.workdir} "
          f"(part {args.part} of {args.nparts})")
    done = fd.run(part=args.part, nparts=args.nparts)
    print(f"computed {done}; status: {fd._count()}")


def cmd_status(args) -> None:
    from .hessian import FiniteDifferenceHessian

    fd = FiniteDifferenceHessian.open(args.workdir)
    if args.clear_locks:
        print(f"removed {fd.clear_locks()} lock(s)")
    st = fd.status()
    print(fd._count())
    if st["running"]:
        print("running/locked:", " ".join(st["running"]))


def cmd_modes(args) -> None:
    from .hessian import FiniteDifferenceHessian

    fd = FiniteDifferenceHessian.open(args.workdir)
    modes = fd.normal_modes()
    f0 = fd.reference_forces()
    print(f"max |force| at the reference geometry: "
          f"{np.abs(f0).max():.2e} eV/Angstrom")
    print(f"max asymmetry of the Hessian: {fd.asymmetry:.2e} eV/Angstrom^2")
    print(modes.summary())
    if len(modes.imaginary):
        print("WARNING: imaginary modes; the reference is not a minimum")
    modes.save(args.modes)
    print(f"saved {args.modes}")


def cmd_lineshape(args) -> None:
    from ase.io import read

    from .coupling import couplings_from_gradient, read_gap_gradient
    from .lineshape import ensemble_spectrum, vibronic_lineshape
    from .modes import EV2INVCM, NormalModes, read_orca_hessian

    if args.modes:
        modes = NormalModes.load(args.modes)
    else:
        if not args.xyz:
            sys.exit("--orca-hess needs --xyz")
        modes = NormalModes.from_hessian(read(args.xyz),
                                         read_orca_hessian(args.orca_hess))
    g = read_gap_gradient(args.gradient)
    c = couplings_from_gradient(modes, g, min_wavenumber=args.min_wavenumber)
    print(c.summary())
    if args.couplings:
        np.savetxt(args.couplings,
                   np.c_[c.wavenumbers, c.hbar_omega, c.reorganization,
                         c.huang_rhys],
                   header="cm^-1 hbar_omega/eV lambda/eV S")
    below = None if args.classical_below is None else \
        args.classical_below / EV2INVCM
    ls = vibronic_lineshape(c, args.temperature, args.process,
                            lambda_classical=args.lambda_classical,
                            sigma=args.sigma, gamma=args.gamma,
                            classical_below=below,
                            step=args.step)
    header = (f"lambda_total {ls.reorganization:.6f} eV, E00 - Evert "
              f"{ls.zero_zero:.6f} eV, T {args.temperature} K, "
              f"lambda_classical {args.lambda_classical} eV, "
              f"sigma {args.sigma} eV, gamma {args.gamma} eV, "
              f"process {args.process}\n")
    if args.sites:
        data = np.loadtxt(args.sites, ndmin=2)
        keep = np.ones(len(data), dtype=bool)
        for w in args.where or []:
            col, op, val = _parse_where(w)
            keep &= op(data[:, col], val)
        sites = data[keep, args.column]
        weights = None if args.weight_column is None else \
            data[keep, args.weight_column]
        e, inten = ensemble_spectrum(sites, ls, weights=weights)
        header += f"ensemble of {keep.sum()} sites\nE/eV intensity"
    else:
        e, inten = ls.energy, ls.intensity
        mask = inten > 1e-10 * inten.max()
        lo, hi = np.flatnonzero(mask)[[0, -1]]
        e, inten = e[lo:hi + 1], inten[lo:hi + 1]
        header += "E-Evert/eV F/eV^-1"
    np.savetxt(args.output, np.c_[e, inten], header=header)
    print(f"lambda_total = {ls.reorganization:.4f} eV; wrote {args.output}")


def run(args) -> None:
    """Run the subcommand of parsed xtp_vibronic arguments."""
    {"hessian": cmd_hessian, "status": cmd_status, "modes": cmd_modes,
     "lineshape": cmd_lineshape}[args.command](args)
