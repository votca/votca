#!/usr/bin/env python3
# """Multi purpose script for HNC inverse methods"""
#
# Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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
# Symbols:
# g: RDF
# U: potential
# dU: potential update (U_{k+1} - U_k)
#
# suffixes:
# _cur: current (of step k if currently doing iteration k)
# _tgt: target
# _ce: core_end (where RDF becomes > 0)
# _co: cut_off
# _sc: single_component
#
# prefixes:
# ndx_: index


import argparse
import sys
import numpy as np


BAR_PER_MD_PRESSURE = 16.6053904
G_MIN = 1e-10
G_MIN_EXTRAPOLATE = 1e-1
OMEGA_SAFE_SCALE = 0.99999
np.seterr(all='raise')


def readin_table(filename):
    """read in votca table"""
    table_dtype = {'names': ('x', 'y', 'y_flag'),
                   'formats': ('f', 'f', 'U2')}
    x, y, y_flag = np.loadtxt(filename, dtype=table_dtype, comments=['#', '@'],
                              unpack=True)
    return x, y, y_flag


def saveto_table(filename, x, y, y_flag, comment=""):
    """save votca table"""
    data = np.zeros((len(x),), dtype='f, f, U2')
    data['f0'] = x
    data['f1'] = y
    data['f2'] = y_flag
    np.savetxt(filename, data, header=comment, fmt=['%e', '%e', '%s'])


def compare_grids(grid_a, grid_b):
    """check two grids for no point differing more than 0.1 picometer"""
    if np.any(grid_a - grid_b > 0.0001):
        raise Exception("Different grids!")


def calc_grid_spacing(grid, relative_tolerance=0.01):
    """returns the spacing of an equidistant 1D grid.
    fails if not equidistant grid."""
    diffs = np.diff(grid)
    if abs((max(diffs) - min(diffs)) / max(diffs)) > relative_tolerance:
        raise Exception('the grid is not equidistant')
    return np.mean(diffs)


def fourier(r, y, omega):
    """fourier transform real data y on grid r into reciprocal space omega"""
    Delta_r = calc_grid_spacing(r)
    y_hat = np.zeros_like(omega)
    with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
        for i, omega_i in enumerate(omega):
            y_hat[i] = (2 / omega_i * Delta_r
                        * np.sum(r * y * np.sin(2 * np.pi * omega_i * r)))
    return y_hat


def gen_omega(r):
    """reciprocal space grid from real space grid"""
    Delta_r = calc_grid_spacing(r)
    nyquist_freq = 1 / Delta_r / 2
    omega = np.linspace(nyquist_freq/len(r), nyquist_freq, len(r))
    omega *= OMEGA_SAFE_SCALE
    return omega


def gen_fourier_matrix(r, omega, fourier_function):
    """make a fourier matrix"""
    fourier_matrix = np.identity(len(omega))
    for col_index, col in enumerate(fourier_matrix.T):
        fourier_matrix.T[col_index] = fourier_function(r, col, omega)
    return fourier_matrix


def find_nearest_ndx(array, value):
    """find index of array where closest to value"""
    array = np.asarray(array)
    ndx = (np.abs(array - value)).argmin()
    return ndx


def find_after_cut_off_ndx(array, cut_off):
    """Find index of array after given cut_off. Assumes array is sorted.
    Used for finding first index after cut_off."""
    array = np.asarray(array)
    ndx_closest = find_nearest_ndx(array, cut_off)
    if np.isclose(array[ndx_closest], cut_off):
        return ndx_closest + 1
    elif array[-1] < cut_off:
        return len(array)
    else:
        ndx = np.where(array > cut_off)[0][0]
        return ndx


def calc_slices(r, g_tgt, g_cur, cut_off, verbose=False):
    # there are different regions in r used
    #              |       crucial     |                     # regions (slices)
    # |   core     |                 nocore               |
    # Δr--------core_end------------cut_off-----------r[-1]  # distances
    # 0----------ndx_ce--------------ndx_co----------len(r)  # indices
    #
    # nocore equals Δ from Delbary et al.
    # crucial equals Δ' from Delbary et al.
    # note: Vector w of HNCGN is in region crucial, but with one element less, because
    # of the antiderivative operator A0
    r, g_tgt, g_cur = map(np.asarray, (r, g_tgt, g_cur))
    ndx_ce = max((np.where(g_tgt > G_MIN)[0][0],
                  np.where(g_cur > G_MIN)[0][0]))
    ndx_co = find_after_cut_off_ndx(r, cut_off)
    crucial = slice(ndx_ce, ndx_co)
    nocore = slice(ndx_ce, len(r))
    if verbose:
        print(f"ndx_ce: {ndx_ce}, ({r[ndx_ce]})")
        print(f"ndx_co: {ndx_co}, ({cut_off})")
        print(f"min(r): {min(r)}")
        print(f"max(r): {max(r)}")
        print(f"len(r): {len(r)}")
        print("crucial:", crucial.start, crucial.stop, min(r[crucial]), max(r[crucial]))
        print("nocore:", nocore.start, nocore.stop, min(r[nocore]), max(r[nocore]))
    return nocore, crucial


def calc_U(r, g_tgt, G_minus_g, n, kBT, rho, closure):
    """Calculate a potential U using integral equation theory. Supports symmetric
    molecules with n equal beads.

    Args:
        r: Distance grid.
        g_tgt: Target RDF.
        G_minus_g: Target intramolecular RDF. Can be an empty array if n == 1.
        n: Number of equal beads per molecule.
        kBT: Boltzmann constant times temperature.
        rho: Number density of the molecules.
        closure: OZ-equation closure ('hnc' or 'py').

    Returns:
        The calculated potential.
    """
    # remove r = 0 for numerical stability
    r0_removed = False
    if np.isclose(r[0], 0.0):
        r, g_tgt, G_minus_g = r[1:], g_tgt[1:], G_minus_g[1:]
        r0_removed = True
    # reciprocal space ω with radial symmetry
    omega = gen_omega(r)
    # total correlation function h
    h = g_tgt - 1
    # FT of h
    h_hat = fourier(r, h, omega)
    # direct correlation function c from OZ
    if n == 1:
        # single bead case
        c_hat = h_hat / (1 + rho * h_hat)
    else:
        G_minus_g_hat = fourier(r, G_minus_g, omega)
        c_hat = h_hat / ((1 + n * rho * G_minus_g_hat)**2
                         + (1 + n * rho * G_minus_g_hat) * n * rho * h_hat)
    c = fourier(omega, c_hat, r)
    # U from HNC
    with np.errstate(divide='ignore', invalid='ignore'):
        if closure == 'hnc':
            U = kBT * (-np.log(g_tgt) + h - c)
        elif closure == 'py':
            U = kBT * np.log(1 - c/g_tgt)
    if r0_removed:
        U = np.insert(U, 0, np.nan)
    return U


def calc_dU_newton(r, g_tgt, g_cur, G_minus_g, n, kBT, rho, cut_off,
                   closure, newton_mod, verbose=False):
    """Calculate a potential update dU using Newtons method and integral equation
    theory. Supports symmetric molecules with n equal beads.

    Args:
        r: Distance grid.
        g_tgt: Target RDF.
        g_cur: Current RDF.
        G_minus_g: Current intramolecular RDF. Can be an empty array if n == 1.
        n: Number of equal beads per molecule.
        kBT: Boltzmann constant times temperature.
        rho: Number density of the molecules.
        cut_off: Highest distance for potential update.
        closure: OZ-equation closure ('hnc' or 'py').
        newton_mod: Use IBI style update term.

    Returns:
        The calculated potential update.
    """
    # remove r = 0 for numerical stability
    r0_removed = False
    if np.isclose(r[0], 0.0):
        r, g_cur, g_tgt, G_minus_g = r[1:], g_cur[1:], g_tgt[1:], G_minus_g[1:]
        r0_removed = True
    # calc slices and Delta_r
    nocore, crucial = calc_slices(r, g_tgt, g_cur, cut_off, verbose=verbose)
    # reciprocal space with radial symmetry ω
    omega = gen_omega(r)
    # difference of rdf to target
    Delta_g = g_cur - g_tgt
    # FT of total correlation function 'h'
    h_hat = fourier(r, g_cur - 1, omega)
    F = gen_fourier_matrix(r, omega, fourier)
    # dc/dg
    if n == 1:
        # single bead case
        dcdg = np.matmul(np.matmul(np.linalg.inv(F),
                                   np.diag(1 / (1 + rho * h_hat)**2)),
                         F)
    else:
        G_minus_g_hat = fourier(r, G_minus_g, omega)
        dcdg = np.matmul(np.matmul(np.linalg.inv(F),
                                   np.diag(1 / (1 + n * rho * G_minus_g_hat
                                                + n * rho * h_hat)**2)),
                         F)
    # calculate jacobian^-1
    if closure == 'hnc':
        if newton_mod:
            with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
                jac_inv = kBT * (np.diag(1 - np.log(g_tgt / g_cur) / Delta_g) + dcdg)
        else:
            with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
                jac_inv = kBT * (np.diag(1 - 1 / g_cur) - dcdg)
    elif closure == 'py':
        raise NotImplementedError
    # cut jacobian and transform back
    jac = np.linalg.inv(jac_inv)
    jac_cut = jac[crucial, crucial]
    jac_inv_cut = np.linalg.inv(jac_cut)
    dU = - np.matmul(jac_inv_cut, Delta_g[crucial])
    if verbose:
        np.savez_compressed('newton-arrays.npz', jac=jac, jac_cut=jac_cut,
                            jac_inv=jac_inv, jac_inv_cut=jac_inv_cut)
    # fill core and behind cut-off
    dU = np.concatenate((np.full(nocore.start + (1 if r0_removed else 0), np.nan),
                         dU,
                         np.full(len(r) - crucial.stop, np.nan)))
    return dU


def gauss_newton_constrained(A, C, b, d):
    """do a gauss-newton update, but eliminate Cx=d first"""
    m, n = A.shape
    p, n = C.shape
    b.shape = (m)
    d.shape = (p)

    if p > 1:
        raise Exception("not implemented for p > 1")

    A_elim = A.copy()
    b_elim = b.copy()
    for i in range(p):
        pivot = np.argmax(abs(C[i]))  # find max value of C
        A_elim = A - (np.ones_like(A) * A[:, pivot][:, np.newaxis]
                      * C[i] / C[i, pivot])
        b_elim = b - A[:, pivot] * d[i] / C[i, pivot]
        A_elim = np.delete(A_elim, pivot, 1)
    if p == n:
        print("WARNING: solution of Gauss-Newton update determined fully "
              "by constraints.")
        x_elim = []
    else:
        x_elim = np.linalg.solve(np.matmul(A_elim.T, A_elim),
                                 np.matmul(A_elim.T, b_elim))
    if p == 0:
        # no constraints
        x = x_elim
    else:
        x_pivot = (d[i] - np.matmul(np.delete(C, pivot, 1),
                                    x_elim)) / C[i, pivot]
        x = np.insert(x_elim, pivot, x_pivot)
    return x


def calc_dU_gauss_newton(r, g_tgt, g_cur, G_minus_g, n, kBT, rho,
                         cut_off, constraints,
                         verbose=False):
    """Calculate a potential update dU using the Gauss-Newton method and
    integral equation theory. Constraints can be added.

    Args:
        r: Distance grid.
        g_tgt: Target RDF.
        g_cur: Current RDF.
        kBT: Boltzmann constant times temperature.
        rho: Number density of the molecules.
        cut_off: Highest distance for potential update.
        constraints: List of dicts, which describe physical constraints.

    Returns:
        The calculated potential update.
    """

    # remove r = 0 for numerical stability
    r0_removed = False
    if np.isclose(r[0], 0.0):
        r, g_cur, g_tgt, G_minus_g = r[1:], g_cur[1:], g_tgt[1:], G_minus_g[1:]
        r0_removed = True

    # calc slices and Delta_r
    nocore, crucial = calc_slices(r, g_tgt, g_cur, cut_off, verbose=verbose)
    Delta_r = calc_grid_spacing(r)
    # reciprocal grid up to Nyquist frequency
    omega = gen_omega(r)
    # Fourier matrix
    F = gen_fourier_matrix(r, omega, fourier)
    # pair correlation function 'h'
    h = g_cur - 1
    # special Fourier of h
    h_hat = fourier(r, h, omega)
    # dc/dg
    if n == 1:
        # single bead case
        dcdg = np.matmul(np.matmul(np.linalg.inv(F),
                                   np.diag(1 / (1 + rho * h_hat)**2)),
                         F)
    else:
        G_minus_g_hat = fourier(r, G_minus_g, omega)
        dcdg = np.matmul(np.matmul(np.linalg.inv(F),
                                   np.diag(1 / (1 + n * rho * G_minus_g_hat
                                                + n * rho * h_hat)**2)),
                         F)
    # jacobian^-1 (matrix U in Delbary et al., with respect to potential)
    with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
        jac_inv = kBT * (np.diag(1 - 1 / g_cur[nocore]) - dcdg[nocore, nocore])
    # A0 matrix
    A0 = Delta_r * np.triu(np.ones((len(r[nocore]), len(r[crucial])-1)), k=0)
    # Jacobian with respect to force
    J = np.matmul(np.linalg.inv(jac_inv), A0)
    # constraint matrix and vector
    C = np.zeros((len(constraints), len(r[crucial])-1))
    d = np.zeros(len(constraints))
    # build constraint matrix and vector from constraints
    if verbose:
        print(constraints)
    for c, constraint in enumerate(constraints):
        if constraint['type'] == 'pressure':
            # current pressure
            p = constraint['current'] / BAR_PER_MD_PRESSURE
            # target pressure
            p_tgt = constraint['target'] / BAR_PER_MD_PRESSURE
            # g_tgt(r_{i+1})
            g_tgt_ip1 = g_tgt[crucial][1:]
            # g_tgt(r_{i})
            g_tgt_i = g_tgt[crucial][:-1]
            # r_{i+1}
            r_ip1 = r[crucial][1:]
            # r_{i}
            r_i = r[crucial][:-1]
            # l vector
            ll = (g_tgt_i + g_tgt_ip1) * (r_ip1**4 - r_i**4)
            ll *= 1/12 * np.pi * rho**2
            # set C row and d element
            C[c, :] = ll
            d[c] = p_tgt - p
        else:
            raise Exception("not implemented constraint type")
    # residuum vector
    res = g_tgt - g_cur
    # switching to notation of Gander et al. for solving
    A = J
    b = res[nocore]
    w = gauss_newton_constrained(A, C, b, d)
    # dU
    dU = np.matmul(A0, w)
    # fill core with nans
    dU = np.concatenate((np.full(nocore.start + (1 if r0_removed else 0), np.nan), dU))
    # dump files
    if verbose:
        np.savez_compressed('gauss-newton-arrays.npz', A=A, b=b, C=C, d=d,
                            jac_inv=jac_inv, A0=A0, J=J)
    return dU


def extrapolate_U_constant(dU, dU_flag):
    """Extrapolate the potential in the core region by a constant value.

    Returns dU. U_{k+1} = U_k + dU is done by Votca."""
    dU_extrap = dU.copy()
    # find first valid dU value
    first_dU_index = np.where(dU_flag == 'i')[0][0]
    first_dU = dU[first_dU_index]

    # replace out of range dU values with constant first value
    dU_extrap = np.where(dU_flag == 'i', dU, first_dU)
    return dU_extrap


def extrapolate_U_power(r, dU, U, g_tgt, g_min, kBT, verbose=False):
    """Extrapolate the potential in the core region including
    the point, where the RDF becomes larger than g_min or where
    the new potential is convex. A power function is used.
    The PMF is fitted, not U+dU. The fit is then shifted such
    that the graph is monotonous.

    Extrapolation is done, because p-HNCGN often has artifacs at
    the core, especially when pressure is far off.

    Returns dU. U_{k+1} = U_k + dU is done by Votca."""
    # make copy
    dU_extrap = dU.copy()
    # calc PMF
    with np.errstate(divide='ignore', over='ignore'):
        pmf = np.nan_to_num(-kBT * np.log(g_tgt))
    # index first minimum
    ndx_fm = np.where(np.nan_to_num(np.diff(pmf)) > 0)[0][0]
    # index core end
    ndx_ce = np.where(g_tgt > g_min)[0][0]
    if verbose:
        print('ndx_fm, r_fm', ndx_fm, r[ndx_fm])
        print('ndx_ce, r_ce', ndx_ce, r[ndx_ce])
    # fit pmf region
    fit_region = slice(ndx_ce, ndx_ce + 3)
    # fit pmf with power function a*x^b
    pmf_shift = -pmf[ndx_fm] + 0.01
    fit_x = np.log(r[fit_region])
    fit_y = np.log(pmf[fit_region] + pmf_shift)
    b, log_a = np.polyfit(fit_x, fit_y, 1)
    a = np.exp(log_a)
    if verbose:
        print('pmf fit a*x^b. Coefficients a, b:', a, b)
    with np.errstate(divide='ignore', over='ignore'):
        pmf_fit = np.nan_to_num(a * r**b - pmf_shift)

    # region to extrapolate
    ndx_ex1 = ndx_ce + 1
    ndx_ex2 = np.where(np.nan_to_num(np.diff(np.diff(U + dU))) > 0)[0][0]
    ndx_ex = max(ndx_ex1, ndx_ex2)
    if verbose:
        print('extrapolate up to:', r[ndx_ex])
    # extrapolate
    U_extrap = U + dU
    U_extrap[:ndx_ex] = pmf_fit[:ndx_ex] + (U_extrap[ndx_ex] - pmf_fit[ndx_ex])
    dU_extrap = U_extrap - U
    return dU_extrap


def shift_U_cutoff_zero(dU, r, U, cut_off):
    """Make potential zero at and beyond cut-off"""
    dU_shift = dU.copy()
    # shift dU to be zero at cut_off and beyond
    ndx_co = find_after_cut_off_ndx(r, cut_off)
    U_before_cut_off = U[ndx_co-1] + dU[ndx_co-1]
    dU_shift -= U_before_cut_off
    dU_shift[ndx_co:] = -1 * U[ndx_co:]
    return dU_shift


def fix_U_near_cut_off_full(r, U, cut_off):
    """Modify the potential close to the cut-off in
    a way, such that it is more smooth. The derivative
    of the potential between the last two points will
    be equal to the derivative between the two points
    before. The original last two points of dU are
    therefore ignored.

    This also helps agains an artifact of p-HNCGN,
    where the last value of dU is a spike."""
    U_fixed = U.copy()
    ndx_co = find_after_cut_off_ndx(r, cut_off)
    second_last_deriv = U[ndx_co-2] - U[ndx_co-3]
    shift = -1.0 * second_last_deriv - U[ndx_co-2]
    # modify up to second last value
    U_fixed[:ndx_co-1] += shift
    return U_fixed


def upd_flag_g_smaller_g_min(flag, g, g_min):
    """Take a flag list, copy it, and set the flag to 'o'utside if g is smaller
    g_min"""
    flag_new = flag.copy()
    for i, gg in enumerate(g):
        if gg < g_min:
            flag_new[i] = 'o'
    return flag_new


def upd_flag_by_other_flag(flag, other_flag):
    """Take a flag list, copy it, and set the flag to 'o'utside where some
    other flag list is 'o'"""
    flag_new = flag.copy()
    for i, of in enumerate(other_flag):
        if of == 'o':
            flag_new[i] = 'o'
    return flag_new


def upd_U_const_first_flag_i(U, flag):
    """Take a potential list, copy it, and set the potential at places where
    flag is 'o' to the first value where flag is 'i'"""
    U_new = U.copy()
    # find first valid U value
    first_U_index = np.where(flag == 'i')[0][0]
    first_U = U_new[first_U_index]
    # replace out of range U values
    U_new = np.where(flag == 'i', U_new, first_U)
    return U_new


def upd_U_zero_beyond_cut_off(r, U, cut_off):
    """Take a potential list, copy it, and shift U to be zero at cut_off and
    beyond"""
    U_new = U.copy()
    index_cut_off = np.searchsorted(r, cut_off)
    U_cut_off = U_new[index_cut_off]
    U_new -= U_cut_off
    U_new[index_cut_off:] = 0
    return U_new


def main():
    # always raise error

    description = """\
    This script calculatess U or ΔU with the HNC methods.
    """
    parser = argparse.ArgumentParser(description=description)
    # subparsers
    subparsers = parser.add_subparsers(dest='subcommand')
    parser_pot_guess = subparsers.add_parser(
        'potential_guess',
        help='potential guess from inverting integral equation')
    parser_newton = subparsers.add_parser(
        'newton',
        help='potential update using Newton method')
    parser_newton_mod = subparsers.add_parser(
        'newton-mod',
        help='potential update using a modified Newton method')
    parser_gauss_newton = subparsers.add_parser(
        'gauss-newton',
        help='potential update using Gauss-Newton method')

    # all subparsers
    for pars in [parser_pot_guess, parser_newton, parser_newton_mod,
                 parser_gauss_newton]:
        pars.add_argument('-v', '--verbose', dest='verbose',
                          help='save some intermeditary results',
                          action='store_const', const=True, default=False)
        pars.add_argument('--closure', type=str, choices=['hnc', 'py'],
                          required=True,
                          help='Closure equation to use for the OZ equation')
        pars.add_argument('--g-tgt', type=argparse.FileType('r'),
                          nargs='+', required=True,
                          metavar=('X-X.dist.tgt', 'X-Y.dist.tgt'),
                          help='RDF target files')
        pars.add_argument('--kBT', type=float, required=True, help='')
        pars.add_argument('--densities', type=float,
                          nargs='+', required=True,
                          metavar=('rho_X', 'rho_Y'),
                          help='list of number densities of beads')
        pars.add_argument('--n-intra', type=int,
                          nargs='+', required=True,
                          metavar=('n_X', 'n_Y'),
                          help='number of beads per molecule')
        # todo: this might need to be split up for multicomponent
        pars.add_argument('--cut-off', type=float, required=True,
                          help='cutt-off (co) of potential')
        pars.add_argument('--U-out', type=argparse.FileType('wb'),
                          nargs='+', required=True,
                          metavar=('X-X.dpot.new', 'X-Y.dpot.new'),
                          help='U or ΔU output files')
    # intial potential guess subparsers
    for pars in [parser_pot_guess]:
        pars.add_argument('--G-tgt', type=argparse.FileType('r'),
                          nargs='+', required=False,
                          metavar=('X-X-incl.dist.tgt', 'X-Y-incl.dist.tgt'),
                          help='RDF (including intramolecular) target files')
    # update potential subparsers
    for pars in [parser_newton, parser_newton_mod, parser_gauss_newton]:
        pars.add_argument('--g-cur', type=argparse.FileType('r'),
                          nargs='+', required=True,
                          metavar=('X-X.dist.cur', 'X-Y.dist.cur'),
                          help='RDF current files')
        pars.add_argument('--G-cur', type=argparse.FileType('r'),
                          nargs='+', required=False,
                          metavar=('X-X-incl.dist.cur', 'X-Y-incl.dist.cur'),
                          help='RDF (including intramolecular) current files')
        pars.add_argument('--U-cur', type=argparse.FileType('r'),
                          nargs='+', required=True,
                          metavar=('X-X.pot.cur', 'X-Y.pot.cur'),
                          help='potential current files')
    # HNCGN only options
    parser_gauss_newton.add_argument('--pressure-constraint',
                                     dest='pressure_constraint',
                                     type=str, default=None)
    parser_gauss_newton.add_argument('--extrap-near-core',
                                     dest='extrap_near_core',
                                     type=str, choices=['none', 'constant',
                                                        'power'])
    parser_gauss_newton.add_argument('--fix-near-cut-off',
                                     dest='fix_near_cut_off',
                                     type=str, choices=['none', 'full-deriv'])

    args = parser.parse_args()

    # check for subcommand
    if args.subcommand is None:
        parser.print_help()
        raise Exception("subcommand needed")

    # close writable files directly due to weird bug, where np.savetxt would
    # write empty file, use filename later.
    # I should report it, but on the other side, argparse opened files do not
    # get closed at any point, so this is better
    for f in args.U_out:
        f.close()
    args.U_out = [f.name for f in args.U_out]

    # infering variables from input
    n_beads = len(args.densities)
    # nr. of elements in triangular matrix incl. diagonal
    n_interactions = (n_beads * (n_beads + 1)) // 2

    # some checks on input
    # test for same number of interactions
    file_arguments_names = ['g_tgt', 'G_tgt', 'g_cur', 'G_cur', 'U_cur',
                            'U_out']
    file_arguments = [(argname, vars(args)[argname])
                      for argname in file_arguments_names
                      if vars(args).get(argname) is not None]
    file_arguments_wrong = [(argname, flist)
                            for argname, flist in file_arguments
                            if len(flist) != n_interactions]
    for argname, flist in file_arguments_wrong:
        raise Exception("""N = {} densities provided, therefore
                        there should be (N * (N + 1)) // 2 = {}
                        files for {}, but {} was
                        provided""".format(n_beads, n_interactions, argname,
                                           [f.name for f in flist]))
    # multicomponent not implemented
    if any((len(files) != 1 for files in [args.g_tgt, args.densities])):
        raise Exception('not implemented for multiple components!')
    if len(args.n_intra) > 1:
        raise Exception('not implemented for multiple components!')
    # if n_intra > 1, also needs G
    # NOT TRUE: may just want to initiate without G
    """
    G_arguments_names = ['G_tgt', 'G_cur']
    if args.n_intra[0] > 1 and not any(vars(args).get(G_name)
                                       for G_name in G_arguments_names):
        raise Exception('If n_intra is larger than 1, you should also '
                        'provide some G')
    """
    # todo: if n_intra == 1, check if G close to g at G[-1]
    # todo for multicomponent: check order of input and output by filename
    # todo for multicomponent: allow not existing X-Y? particles would overlap
    # todo for multicomponent: do not allow same bead on different moleculetypes

    # load input arrays
    input_arrays = {}  # input_arrays['g_tgt'][0]['x']
    input_arguments_names = ['g_tgt', 'G_tgt', 'g_cur', 'G_cur', 'U_cur']
    input_arguments = [(argname, vars(args)[argname])
                       for argname in input_arguments_names
                       if vars(args).get(argname) is not None]
    for argname, flist in input_arguments:
        input_arrays[argname] = []
        for i, f in enumerate(flist):
            x, y, flag = readin_table(f)
            input_arrays[argname].append({'x': x, 'y': y, 'flag': flag})

    # todo: compare grids of all
    r = input_arrays['g_tgt'][0]['x']
    # compare_grids(r, input_arrays['g_cur'][0]['x'])
    # compare grids of G_cur and G_tgt with g_tgt in smart way

    # calculate G_minus_g
    if 'G_tgt' in input_arrays:
        g_name = 'g_tgt'
        G_name = 'G_tgt'
    elif 'G_cur' in input_arrays:
        g_name = 'g_cur'
        G_name = 'G_cur'
    else:
        print('no intramolecular correlations provided, assuming there are '
              'none')
        # this will result in zero arrays in input_arrays['G_minus_g'] in the
        # below loop
        g_name = 'g_tgt'
        G_name = 'g_tgt'
    input_arrays['G_minus_g'] = []
    for g_dict, G_dict in zip(input_arrays[g_name], input_arrays[G_name]):
        G_minus_g = np.zeros_like(g_dict['y'])
        G_minus_g[:G_dict['y'].shape[0]] = (G_dict['y']
                                            - g_dict['y']
                                            [:G_dict['y'].shape[0]])
        input_arrays['G_minus_g'].append({'y': G_minus_g})

    # guess potential from distribution
    if args.subcommand in ['potential_guess']:
        U1 = calc_U(r,
                    input_arrays['g_tgt'][0]['y'],
                    input_arrays['G_minus_g'][0]['y'],
                    args.n_intra[0], args.kBT, args.densities[0],
                    args.closure)
        U_flag1 = np.array(['i'] * len(r))
        U_flag2 = upd_flag_g_smaller_g_min(U_flag1,
                                           input_arrays['g_tgt'][0]['y'],
                                           G_MIN)
        U2 = upd_U_const_first_flag_i(U1, U_flag2)
        U3 = upd_U_zero_beyond_cut_off(r, U2, args.cut_off)
        comment = "created by: {}".format(" ".join(sys.argv))
        saveto_table(args.U_out[0], r, U3, U_flag2, comment)

    # newton update
    if args.subcommand in ['newton', 'newton-mod']:
        dU1 = calc_dU_newton(r,
                             input_arrays['g_tgt'][0]['y'],
                             input_arrays['g_cur'][0]['y'],
                             input_arrays['G_minus_g'][0]['y'],
                             args.n_intra[0], args.kBT,
                             args.densities[0], args.cut_off, args.closure,
                             args.subcommand == 'newton-mod',
                             verbose=args.verbose)
        dU_flag1 = np.array(['i'] * len(r))
        dU_flag2 = upd_flag_g_smaller_g_min(dU_flag1,
                                            input_arrays['g_tgt'][0]['y'],
                                            G_MIN)
        dU_flag3 = upd_flag_g_smaller_g_min(dU_flag2,
                                            input_arrays['g_cur'][0]['y'],
                                            G_MIN)
        dU_flag4 = upd_flag_by_other_flag(dU_flag3,
                                          input_arrays['U_cur'][0]['flag'])

        dU2 = upd_U_const_first_flag_i(dU1, dU_flag4)
        U_temp = upd_U_zero_beyond_cut_off(r,
                                           dU2 + input_arrays['U_cur'][0]['y'],
                                           args.cut_off)
        dU3 = U_temp - input_arrays['U_cur'][0]['y']
        comment = "created by: {}".format(" ".join(sys.argv))
        saveto_table(args.U_out[0], r, dU3, dU_flag4, comment)

    # gauss-newton update
    if args.subcommand in ['gauss-newton']:
        # parse constraints
        constraints = []
        if args.pressure_constraint is not None:
            p_target = float(args.pressure_constraint.split(',')[0])
            p_current = float(args.pressure_constraint.split(',')[1])
            constraints.append({'type': 'pressure', 'target': p_target,
                                'current': p_current})

        # calc dU_pure
        dU_pure = calc_dU_gauss_newton(r,
                                       input_arrays['g_tgt'][0]['y'],
                                       input_arrays['g_cur'][0]['y'],
                                       input_arrays['G_minus_g'][0]['y'],
                                       args.n_intra[0], args.kBT,
                                       args.densities[0], args.cut_off,
                                       constraints, verbose=args.verbose)

        # set dU_flag to 'o' inside the core
        dU_flag = np.where(np.isnan(dU_pure), 'o', 'i')

        # select extrapolation
        if args.extrap_near_core == 'none':
            dU_extrap = np.nan_to_num(dU_pure)
        elif args.extrap_near_core == 'constant':
            dU_extrap = extrapolate_U_constant(dU_pure, dU_flag)
        elif args.extrap_near_core == 'power':
            dU_extrap = extrapolate_U_power(r, dU_pure,
                                            input_arrays['U_cur'][0]['y'],
                                            input_arrays['g_tgt'][0]['y'],
                                            G_MIN_EXTRAPOLATE, args.kBT,
                                            verbose=args.verbose)
        else:
            raise Exception("unknown extrapolation scheme for inside and near "
                            "core region: " + args.extrap_near_core)
        # shifts to correct potential after cut-off
        dU_shift = shift_U_cutoff_zero(dU_extrap, r,
                                       input_arrays['U_cur'][0]['y'],
                                       args.cut_off)
        # shifts to correct potential near cut-off
        if args.fix_near_cut_off == 'none':
            dU = dU_shift.copy()
        elif args.fix_near_cut_off == 'full-deriv':
            U_new = input_arrays['U_cur'][0]['y'] + dU_shift
            U_new = fix_U_near_cut_off_full(r, U_new, args.cut_off)
            dU = U_new - input_arrays['U_cur'][0]['y']
        else:
            raise Exception("unknown fix scheme for near cut-off: "
                            + args.fix_near_cut_off)

        if args.verbose:
            np.savez_compressed('hncgn-dU.npz', r=r, dU_pure=dU_pure,
                                dU_extrap=dU_extrap, dU_shift=dU_shift)
        comment = "created by: {}".format(" ".join(sys.argv))
        saveto_table(args.U_out[0], r, dU, dU_flag, comment)


if __name__ == '__main__':
    main()
