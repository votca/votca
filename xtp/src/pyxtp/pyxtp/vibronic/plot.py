"""Publication-quality figures of spectra and vibronic couplings.

Needs matplotlib (not imported by ``pyxtp.vibronic`` itself).

    from pyxtp.vibronic.plot import paper_style, plot_spectrum

    with paper_style():
        fig = plot_spectrum(spectrum, broadening=0.3, onset=True,
                            experiment=(e_exp, i_exp))
        fig.savefig("ipes.pdf")

Design: one axis per panel; thin lines; the contributions of the levels
in a fixed colour order (colour-blind safe; HOMO and LUMO first, then by
distance from the gap, so a level keeps its colour; from the ninth level
on grey), the totals and all text in black/grey ink; a legend
whenever there is more than one curve.
"""

from __future__ import annotations

import contextlib
from typing import Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from .coupling import VibronicCouplings
from .spectrum import Spectrum

#: categorical colours in fixed order (validated for colour-vision
#: deficiencies on a white surface); never cycled
SERIES_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                 "#008300", "#6250d6", "#e34948")
INK = "#1a1a19"
INK_SECONDARY = "#6b6a63"
INK_MUTED = "#a3a299"

#: figure widths in inches: one and two journal columns
SINGLE_COLUMN = 3.375
DOUBLE_COLUMN = 7.0


@contextlib.contextmanager
def paper_style(fontsize: float = 8.0):
    """Matplotlib settings for journal figures (vector fonts, thin lines)."""
    rc = {
        "font.size": fontsize, "axes.labelsize": fontsize,
        "axes.titlesize": fontsize, "legend.fontsize": fontsize - 1,
        "xtick.labelsize": fontsize - 1, "ytick.labelsize": fontsize - 1,
        "font.family": "sans-serif", "mathtext.fontset": "dejavusans",
        "axes.linewidth": 0.6, "axes.edgecolor": INK_SECONDARY,
        "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK_SECONDARY, "ytick.color": INK_SECONDARY,
        "xtick.labelcolor": INK, "ytick.labelcolor": INK,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.0, "lines.solid_capstyle": "round",
        "legend.frameon": False, "legend.handlelength": 1.6,
        "figure.dpi": 150, "savefig.dpi": 600,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    }
    with plt.rc_context(rc):
        yield


def _slot(label: str, index: int) -> int:
    """Colour slot of a contribution: n for HOMO-n and LUMO+n."""
    for name, sep in (("HOMO", "-"), ("LUMO", "+")):
        if label == name:
            return 0
        if label.startswith(name + sep) and label[5:].isdigit():
            return int(label[5:])
    return index


def _new_axes(ax, width, height):
    if ax is not None:
        return ax.figure, ax
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    return fig, ax


def plot_spectrum(spectrum: Spectrum,
                  broadening: Optional[float] = None,
                  experiment: Optional[Tuple[Sequence[float],
                                             Sequence[float]]] = None,
                  onset: bool = False,
                  contributions: bool = True,
                  xlim: Optional[Tuple[float, float]] = None,
                  experiment_label: str = "Experiment",
                  ax=None, width: float = SINGLE_COLUMN,
                  height: Optional[float] = None):
    """Spectrum with the contributions of the levels.

    Args:
        spectrum: Spectrum (e.g. from level_spectrum)
        broadening: standard deviation (eV) of a Gaussian instrument
            broadening; the broadened total is drawn in addition to the
            intrinsic one
        experiment: (energy, intensity) of a measured spectrum; scaled to
            the maximum of the (broadened) total inside xlim
        onset: mark the tangent onset and the first peak of the
            (broadened) total; the values are in fig.onset
        contributions: draw the individual levels
        xlim: energy range
        ax: existing axes (otherwise a new single-column figure)

    Returns:
        the figure; fig.onset holds the onset dictionary if requested
    """
    fig, ax = _new_axes(ax, width, height or min(0.72 * width, 2.6))
    e = spectrum.energy
    total = spectrum.total
    broad = spectrum.broadened(broadening) if broadening else None
    main_e, main_y = (broad.energy, broad.total) if broad else (e, total)
    if xlim is None:
        mask = main_y > 2e-3 * main_y.max()
        xlim = (main_e[mask][0], main_e[mask][-1])
    in_main = (main_e >= xlim[0]) & (main_e <= xlim[1])
    in_e = (e >= xlim[0]) & (e <= xlim[1])
    # common scale: the intrinsic total has maximum 1 inside xlim
    scale = 1.0 / total[in_e].max()
    labels = list(spectrum.contributions)
    many = len(labels) > 1

    if contributions and many:
        # colour by the distance of the level from the gap (HOMO and LUMO
        # first), so that a level keeps its colour whatever else is shown;
        # the legend lists the levels in this order as well
        order = sorted(range(len(labels)),
                       key=lambda j: (_slot(labels[j], j), j))
        grey_label = "deeper/higher levels"
        for i, label in ((_slot(labels[j], j), labels[j]) for j in order):
            y = spectrum.contributions[label] * scale
            if i < len(SERIES_COLORS):
                ax.plot(e, y, color=SERIES_COLORS[i], lw=1.0, label=label,
                        zorder=3)
            else:
                ax.plot(e, y, color=INK_MUTED, lw=0.8, zorder=3,
                        label=grey_label)
                grey_label = None
    if broad:
        # below the contributions: it coincides with them where the levels
        # do not overlap
        ax.plot(e, total * scale, color=INK_SECONDARY, lw=1.0, zorder=2,
                label="Total, intrinsic" if many else "Intrinsic")
        # the broadened curve is rescaled to the same height for comparison
        bscale = 1.0 / main_y[in_main].max()
        ax.plot(main_e, main_y * bscale, color=INK, lw=1.3, zorder=4,
                label=f"Broadened ($\\sigma$ = {broadening:g} eV)")
        top_y = main_y * bscale
    else:
        ax.plot(e, total * scale, color=INK, lw=1.3,
                zorder=2 if (contributions and many) else 4,
                label="Total" if many else None)
        top_y = total * scale
    if experiment is not None:
        xe = np.asarray(experiment[0], dtype=float)
        ye = np.asarray(experiment[1], dtype=float)
        sel = (xe >= xlim[0]) & (xe <= xlim[1])
        ye = ye / ye[sel].max()
        ax.plot(xe, ye, ls="none", marker="o", ms=2.5, mfc="none",
                mec=INK_SECONDARY, mew=0.6, label=experiment_label)
    fig.onset = None
    if onset:
        res = Spectrum(main_e, {"t": top_y}, spectrum.process).onset()
        fig.onset = res
        x0, y0, slope = res["x0"], res["y0"], res["slope"]
        ypk = float(np.interp(res["peak"], main_e, top_y))
        xt = np.array([res["onset"], x0 + (ypk - y0) / slope])
        ax.plot(xt, y0 + slope * (xt - x0), color=INK, lw=0.7, ls=(0, (4, 2)))
        ax.plot([res["onset"]], [0.0], marker="|", ms=7, color=INK, mew=0.9,
                clip_on=False)
        ax.plot([res["peak"]], [ypk], marker="o", ms=3, color=INK,
                mec="white", mew=0.6)
        left = slope > 0
        # room for the onset label beside the edge
        span = xlim[1] - xlim[0]
        room = span * 0.75 / fig.get_figwidth()  # about 0.75 inch
        if left and res["onset"] - xlim[0] < room:
            xlim = (res["onset"] - room, xlim[1])
        elif not left and xlim[1] - res["onset"] < room:
            xlim = (xlim[0], res["onset"] + room)
        ax.annotate(f"onset\n{res['onset']:.2f} eV", (res["onset"], 0.0),
                    xytext=(-4 if left else 4, 6), textcoords="offset points",
                    ha="right" if left else "left", va="bottom")
        ax.annotate(f"{res['peak']:.2f} eV", (res["peak"], ypk),
                    xytext=(0, 5), textcoords="offset points", ha="center",
                    va="bottom")
    ax.set_xlim(*xlim)
    ax.set_ylim(0.0, 1.15)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(spectrum.axis_label)
    ax.set_ylabel("Intensity (arb. units)")
    handles, _ = ax.get_legend_handles_labels()
    if len(handles) > 1:
        # above the axes, so that the plot keeps the full column width
        ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0),
                  ncol=3 if len(handles) > 4 else 2, borderaxespad=0.2,
                  columnspacing=1.0, handletextpad=0.5)
    return fig


def plot_couplings(couplings: VibronicCouplings, quantity: str = "lambda",
                   label_modes: int = 3, width: float = SINGLE_COLUMN,
                   height: Optional[float] = None, title: str = ""):
    """Mode-resolved couplings (top) and their running sum (bottom).

    Args:
        couplings: VibronicCouplings
        quantity: "lambda" (reorganization energy, meV) or "huang-rhys"
        label_modes: label this many of the strongest modes with their
            frequency
    """
    if quantity not in ("lambda", "huang-rhys"):
        raise ValueError("quantity must be 'lambda' or 'huang-rhys'")
    wn = couplings.wavenumbers
    order = np.argsort(wn)
    wn = wn[order]
    if quantity == "lambda":
        y = 1000.0 * couplings.reorganization[order]
        ylabel, cumlabel = r"$\lambda_\nu$ (meV)", r"$\sum\lambda_\nu$ (meV)"
    else:
        y = couplings.huang_rhys[order]
        ylabel, cumlabel = r"$S_\nu$", r"$\sum S_\nu$"
    fig, (top, bottom) = plt.subplots(
        2, 1, sharex=True, figsize=(width, height or 0.85 * width),
        constrained_layout=True, gridspec_kw={"height_ratios": [3, 2]})
    color = SERIES_COLORS[0]
    top.vlines(wn, 0.0, y, color=color, lw=1.2)
    top.set_ylim(0.0, 1.18 * y.max())
    top.set_ylabel(ylabel)
    for k in np.argsort(y)[::-1][:label_modes]:
        top.annotate(f"{wn[k]:.0f}", (wn[k], y[k]), xytext=(0, 2),
                     textcoords="offset points", ha="center", va="bottom",
                     fontsize=plt.rcParams["font.size"] - 1.5)
    if title:
        top.set_title(title, loc="left")
    cum = np.cumsum(y)
    bottom.step(np.concatenate([[0.0], wn, [1.05 * wn.max()]]),
                np.concatenate([[0.0], cum, [cum[-1]]]), where="post",
                color=color, lw=1.2)
    bottom.set_ylim(0.0, 1.12 * cum[-1])
    bottom.set_ylabel(cumlabel)
    bottom.set_xlabel(r"Wavenumber (cm$^{-1}$)")
    bottom.set_xlim(0.0, 1.05 * wn.max())
    bottom.annotate(f"{cum[-1]:.1f}" if quantity == "lambda"
                    else f"{cum[-1]:.2f}", (1.05 * wn.max(), cum[-1]),
                    xytext=(-2, 2), textcoords="offset points", ha="right",
                    va="bottom")
    fig.align_ylabels([top, bottom])
    return fig
