"""Subcommands of the xtp_vibronic script.

The parser is in the script, so that its help works without pyxtp.

Subcommands:

  hessian    run (or continue) the displaced force calculations; several
             processes or array-job parts can share one workdir
  status     progress of a workdir
  modes      collect the Hessian, print the frequencies, save modes.npz
  lineshape  couplings and line shape of one transition (optionally an
             ensemble spectrum)
  gradients  add QP gradients or excited-state forces to the orb file of
             an earlier dftgwbse run (no DFT, no GW-BSE)
  spectrum   UPS/IPES-type spectrum of all QP levels with gradients
  plot       publication-quality figure of a spectrum or of couplings
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
    """Run (or continue) the displaced force calculations."""
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
    """Print the progress of a Hessian workdir."""
    from .hessian import FiniteDifferenceHessian

    fd = FiniteDifferenceHessian.open(args.workdir)
    if args.clear_locks:
        print(f"removed {fd.clear_locks()} lock(s)")
    st = fd.status()
    print(fd._count())
    if st["running"]:
        print("running/locked:", " ".join(st["running"]))


def cmd_modes(args) -> None:
    """Collect the Hessian, print the frequencies, save the modes."""
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
    """Compute the couplings and line shape of one transition."""
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


def _load_modes(args):
    from ase.io import read

    from .modes import NormalModes, read_orca_hessian

    if args.modes:
        return NormalModes.load(args.modes)
    if not args.xyz:
        sys.exit("--orca-hess needs --xyz")
    return NormalModes.from_hessian(read(args.xyz),
                                    read_orca_hessian(args.orca_hess))


def cmd_spectrum(args) -> None:
    """Compute the spectrum of all QP levels with gradients."""
    from .modes import EV2INVCM
    from .spectrum import QPLevels, level_spectrum

    modes = _load_modes(args)
    levels = QPLevels.from_orb(args.gradients)
    sites = None
    weights = None
    if args.sites:
        if not args.columns:
            sys.exit("--sites needs --columns, e.g. lumo=3,lumo+1=4")
        data = np.loadtxt(args.sites, ndmin=2)
        keep = np.ones(len(data), dtype=bool)
        for w in args.where or []:
            col, op, val = _parse_where(w)
            keep &= op(data[:, col], val)
        sites = {}
        for item in args.columns.split(","):
            if "=" not in item:
                sys.exit(f"--columns expects level=column, got {item}")
            name, col = item.split("=")
            sites[name.strip()] = data[keep, int(col)]
        if args.weight_column is not None:
            weights = data[keep, args.weight_column]
    below = None if args.classical_below is None else \
        args.classical_below / EV2INVCM
    sp = level_spectrum(
        modes, levels, args.process, args.temperature, site_energies=sites,
        weights=weights,
        select=args.levels.split(",") if args.levels else None,
        lambda_classical=args.lambda_classical, sigma=args.sigma,
        gamma=args.gamma, classical_below=below,
        min_wavenumber=args.min_wavenumber, step=args.step)
    print("level      lambda/eV  placed at")
    for label in sp.contributions:
        print(f"{label:9s} {sp.reorganization[label]:9.4f}   "
              f"{sp.notes[label]}")
    sp.save(args.output,
            header=f"T {args.temperature} K, lambda_classical "
                   f"{args.lambda_classical} eV, sigma {args.sigma} eV, "
                   f"gamma {args.gamma} eV")
    print(f"wrote {args.output}")


def cmd_gradients(args) -> None:
    """Add QP gradients or excited-state forces to an orb file."""
    from ..xtp import xtp
    from .spectrum import QPLevels, add_gradients

    if not args.levels and not args.state:
        sys.exit("give --levels and/or --state")
    options = None
    if args.options:
        options = xtp().options
        _apply_options_file(options, args.options)
    add_gradients(args.orb, qp_levels=args.levels, excited_state=args.state,
                  options=options, settings=_settings(args.set),
                  nthreads=args.threads)
    if args.levels:
        print("QP gradients in", args.orb + ":",
              " ".join(QPLevels.from_orb(args.orb).labels))
    if args.state:
        print(f"forces of {args.state} stored in {args.orb}")


def cmd_plot(args) -> None:
    """Plot a spectrum or vibronic couplings."""
    import matplotlib
    matplotlib.use("Agg")

    from .coupling import VibronicCouplings
    from .modes import EV2INVCM
    from .plot import (DOUBLE_COLUMN, SINGLE_COLUMN, paper_style,
                       plot_couplings, plot_spectrum)
    from .spectrum import Spectrum

    width = DOUBLE_COLUMN if args.wide else SINGLE_COLUMN
    with paper_style(args.fontsize):
        if args.couplings:
            d = np.loadtxt(args.couplings, ndmin=2)
            c = VibronicCouplings(d[:, 0] / EV2INVCM, d[:, 2],
                                  np.arange(len(d)))
            fig = plot_couplings(c, quantity=args.quantity, width=width,
                                 title=args.title or "")
        else:
            sp = Spectrum.load(args.spectrum)
            experiment = None
            if args.experiment:
                x = np.loadtxt(args.experiment, ndmin=2)
                experiment = (x[:, 0], x[:, 1])
            fig = plot_spectrum(sp, broadening=args.broadening,
                                experiment=experiment, onset=args.onset,
                                contributions=not args.no_contributions,
                                xlim=tuple(args.xlim) if args.xlim else None,
                                width=width)
            if fig.onset:
                print(f"onset {fig.onset['onset']:.3f} eV, first peak "
                      f"{fig.onset['peak']:.3f} eV")
        fig.savefig(args.output)
    print(f"wrote {args.output}")


def run(args) -> None:
    """Run the subcommand of parsed xtp_vibronic arguments."""
    {"hessian": cmd_hessian, "status": cmd_status, "modes": cmd_modes,
     "lineshape": cmd_lineshape, "spectrum": cmd_spectrum,
     "gradients": cmd_gradients, "plot": cmd_plot}[args.command](args)
