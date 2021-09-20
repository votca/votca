#!/usr/bin/env python3
"""Multi purpose script for Iterative Integral Equation methods."""
#
# Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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
import xml.etree.ElementTree as ET
try:
    import numpy as np
except ImportError:
    print("Numpy is not installed, but needed for the iterative integral equation "
          "methods.")
    raise
if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.5+.")
from csg_functions import (
    readin_table, saveto_table, calc_grid_spacing, fourier, fourier_all,
    gen_beadtype_property_array, gen_fourier_matrix, find_after_cut_off_ndx,
    r0_removal, get_non_bonded, get_density_dict, get_n_intra_dict,
    gen_interaction_matrix, gen_interaction_dict, gauss_newton_constrained,
    upd_flag_g_smaller_g_min, upd_flag_by_other_flag, gen_flag_isfinite,
    extrapolate_dU_left_constant
)


BAR_PER_MD_PRESSURE = 16.6053904
G_MIN = 1e-10
G_MIN_EXTRAPOLATE = 1e-1
np.seterr(all='raise')


def gen_Omega_hat_mat(G_minus_g_hat_mat, rhos, n_intra):
    # σ is any row sum of ω
    sigma_R = G_minus_g_hat_mat @ np.diag(rhos) + np.identity(len(rhos))
    # weighting of row sum σ
    Omega_hat_mat = np.diag(np.sqrt(n_intra)) @ sigma_R @ np.diag(1/np.sqrt(n_intra))
    return Omega_hat_mat


def adapt_reduced_matrix(mat, n_intra):
    """Adapt the prefactors of a matrix to be compatible with the symmetry reduced RISM
    equation.

    The input matrix is already reduced (rows are atom types not atoms), but factors
    need to be applied."""
    Mat = np.diag(np.sqrt(n_intra)) @ mat @ np.diag(np.sqrt(n_intra))
    return Mat


def unadapt_reduced_matrix(Mat, n_intra):
    """Unadapt the prefactors of a matrix compatible with the symmetry reduced RISM
    equation back to the regular form.

    The input matrix is already reduced (rows are atom types not atoms) and adapted.
    Factors are applied to get back to the regular matrix."""
    mat = np.diag(1/np.sqrt(n_intra)) @ Mat @ np.diag(1/np.sqrt(n_intra))
    return mat


def transpose(mat):
    """First dimension is radius or k. Transpose means swapping the last two axis."""
    return np.swapaxes(mat, -1, -2)


def calc_c_matrix(r, k, h_hat_mat, G_minus_g_hat_mat, rhos, n_intra, verbose=False):
    """Calculate the direct correlation function c from g for all interactions."""
    # row sum of ω, after Bertagnolli and my own notes
    Omega_hat_mat = gen_Omega_hat_mat(G_minus_g_hat_mat, rhos, n_intra)
    # H_hat_mat after Bertagnolli
    H_hat_mat = adapt_reduced_matrix(h_hat_mat, n_intra)
    # Rho_mat after Bertagnolli
    Rhos = rhos / n_intra
    # intermediate step
    # have to transpose to solve x·a = b with numpy by solving a'·x' = b'
    H_over_Omega_plus_rho_H = transpose(np.linalg.solve(
        transpose(Omega_hat_mat + np.diag(Rhos) @ H_hat_mat),
        transpose(H_hat_mat)))
    # direct correlation function C from symmetry reduced OZ
    C_hat_mat = np.linalg.solve(Omega_hat_mat, H_over_Omega_plus_rho_H)
    # c_hat from C_hat
    c_hat_mat = unadapt_reduced_matrix(C_hat_mat, n_intra)
    # c from c_hat
    _, c_mat = fourier_all(k, c_hat_mat)
    if verbose:
        np.savez_compressed('calc-c-matrix.npz', r=r, k=k,
                            h_hat_mat=h_hat_mat, G_minus_g_hat_mat=G_minus_g_hat_mat,
                            rhos=rhos, c_mat=c_mat)
    return c_mat


def calc_g(r, c, G_minus_g, n, rho):
    """Calculate the radial distribution function g(r) from c(r)."""
    r0_removed, (r, c, G_minus_g) = r0_removal(r, c, G_minus_g)
    k, c_hat = fourier(r, c)
    if n == 1:
        h_hat = c_hat / (1 - rho * c_hat)
    else:
        _, G_minus_g_hat = fourier(r, G_minus_g)
        h_hat = ((c_hat * (1 + n * rho * G_minus_g_hat)**2)
                 / (1 - n * rho * (1 + n * rho * G_minus_g_hat) * c_hat))
    _, h = fourier(k, h_hat)
    g = h + 1
    if r0_removed:
        g = np.concatenate(([np.nan], g))
    return g


def calc_dc_ext(r_short, r_long, c_k_short, g_k_short, g_tgt_short, G_minus_g_short, n,
                rho):
    """
    Calculate Δc_ext with netwon method.

    This term is used in the iterative extrapolation of g(r). Jacobian has an
    implicit extrapolation of c with zeros on twice its original range.

    """
    # _k is iteration k
    # _s is short
    # _tgt is target
    r0_removed, (r_short, c_k_short, g_k_short, g_tgt_short,
                 G_minus_g_short) = r0_removal(r_short, c_k_short, g_k_short,
                                               g_tgt_short, G_minus_g_short)
    F = gen_fourier_matrix(r_long, fourier)
    Finv = np.linalg.inv(F)
    B = np.concatenate((np.diag(np.ones(len(r_short))),
                        np.zeros((len(r_long) - len(r_short), len(r_short)))),
                       axis=0)
    Binv = np.linalg.pinv(B)
    J = Binv @ (Finv @ np.diag((1 + n * rho * F @ B @ G_minus_g_short)**2
                               / (1 - (1 + n * rho * F @ B @ G_minus_g_short)
                                  * n * rho * F @ B @ c_k_short)**2) @ F) @ B
    Jinv = np.linalg.pinv(J)
    Δc = -1 * Jinv @ (g_k_short - g_tgt_short)
    if r0_removed:
        Δc = np.concatenate(([0], Δc))
    return Δc


def extrapolate_g(r_short, r_long, g_short, G_minus_g_short, n, rho,
                  k_max=5, output_c=False, verbose=False):
    """
    Extrapolate an RDF to larger r with integral equation theory.

    Assumes c = 0 in the extrapolated region. This is not a good aproximation
    in systems with bonds.

    Args:
        r_short: Input grid.
        r_long: Output grid.
        g_short: Input RDF.
        G_minus_g_short: Intramolecular distribution.
        n: Number of equal beads per molecule.
        rho: Number density of the molecules.
        k_max: Number of iterations.
        output_c: Wether to output the final direct correlation function.
        verbose: Print convergence, dump direct correlation function.

    Returns:
        The extrapolated RDF and depending on output_c the c.

    """
    r0_removed, (r_short, r_long, g_short,
                 G_minus_g_short) = r0_removal(r_short, r_long, g_short,
                                               G_minus_g_short)
    ndx_co = len(r_short)
    G_minus_g_long = np.concatenate((G_minus_g_short,
                                     np.zeros(len(r_long) - len(r_short))))
    # starting guess for c
    c_short_k = [calc_c(r_short, g_short, G_minus_g_short, n, rho)]
    c_long_k = [np.zeros_like(r_long)]
    c_long_k[0][:ndx_co] = c_short_k[0]
    # evaluate starting guess
    g_long_k = [calc_g(r_long, c_long_k[0], G_minus_g_long, n, rho)]
    g_short_k = [g_long_k[0][:ndx_co]]
    # Newton iterations
    for it in range(1, k_max+1):
        # update c
        c_short_k.append(c_short_k[-1]
                         + calc_dc_ext(r_short, r_long, c_short_k[-1], g_short_k[-1],
                                       g_short, G_minus_g_short, n, rho))
        c_long_k.append(np.zeros_like(r_long))
        c_long_k[-1][:ndx_co] = c_short_k[-1]
        # new g
        g_long_k.append(calc_g(r_long, c_long_k[-1], G_minus_g_long, n, rho))
        g_short_k.append(g_long_k[-1][:ndx_co])

    if r0_removed:
        for it in range(0, k_max+1):
            c_short_k[it] = np.concatenate(([np.nan], c_short_k[it]))
            c_long_k[it] = np.concatenate(([np.nan], c_long_k[it]))
            g_long_k[it] = np.concatenate(([np.nan], g_long_k[it]))
            g_short_k[it] = np.concatenate(([np.nan], g_short_k[it]))
    if verbose:
        np.savez_compressed('g-extrapolation.npz', c_short_k=c_short_k,
                            c_long_k=c_long_k, g_long_k=g_long_k, g_short_k=g_short_k)
    if output_c:
        return g_long_k[-1], c_long_k[-1]
    return g_long_k[-1]


def calc_slices(r, g_tgt, g_cur, cut_off, verbose=False):
    """
    Generate slices for the regions used in the IIE methods.

    There are different regions used:
                 |       crucial     |                     # regions (slices)
    |   core     |                 nocore               |
    0---------core_end------------cut_off-----------r[-1]  # distances
    0----------ndx_ce--------------ndx_co----------len(r)  # indices
    nocore equals Δ from Delbary et al.
    crucial equals Δ' from Delbary et al.
    note: Vector w of HNCGN is in region crucial, but with one element less,
    because of the antiderivative operator A0.

    """
    r, g_tgt, g_cur = map(np.asarray, (r, g_tgt, g_cur))
    ndx_ce = max((np.where(g_tgt > G_MIN)[0][0],
                  np.where(g_cur > G_MIN)[0][0]))
    ndx_co = find_after_cut_off_ndx(r, cut_off)
    crucial = slice(ndx_ce, ndx_co)
    nocore = slice(ndx_ce, len(r))
    if verbose:
        print("ndx_ce: {}, ({})".format(ndx_ce, r[ndx_ce]))
        print("ndx_co: {}, ({})".format(ndx_co, cut_off))
        print("min(r): {}".format(min(r)))
        print("max(r): {}".format(max(r)))
        print("len(r): {}".format(len(r)))
        print("crucial:", crucial.start, crucial.stop, min(r[crucial]), max(r[crucial]))
        print("nocore:", nocore.start, nocore.stop, min(r[nocore]), max(r[nocore]))
    return nocore, crucial


def calc_U(r, g_tgt, G_minus_g, n, kBT, rho, closure):
    """
    Calculate a potential U using integral equation theory.

    Supports symmetric molecules with n equal beads.

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
    h = g_tgt - 1
    c = calc_c(r, g_tgt, G_minus_g, n, rho)
    with np.errstate(divide='ignore', invalid='ignore'):
        if closure == 'hnc':
            U = kBT * (-np.log(g_tgt) + h - c)
        elif closure == 'py':
            U = kBT * np.log(1 - c/g_tgt)
    return U


def calc_U_matrix(r, k, g_mat, h_hat_mat, G_minus_g_hat_mat, rhos, n_intra, kBT,
                  closure, verbose=False):
    """
    Calculate a potential U using integral equation theory.

    Args:
        r: Distance grid.
        g_mat: matrix of RDF
        h_hat_mat: matrix of Fourier of TCF
        G_minus_g_mat: matrix of Fourier of intramolecular RDF
        rhos: array of densities of the bead types
        n_intra: array with number of bead per molecule
        kBT: Boltzmann constant times temperature.
        closure: OZ-equation closure ('hnc' or 'py').
        verbose: output calc_U_matrix.npz

    Returns:
        matrix of the calculated potentias.
    """
    # calculate direct correlation function
    c_mat = calc_c_matrix(r, k, h_hat_mat, G_minus_g_hat_mat, rhos, n_intra, verbose)
    with np.errstate(divide='ignore', invalid='ignore'):
        if closure == 'hnc':
            U_mat = kBT * (-np.log(g_mat) + (g_mat - 1) - c_mat)
        elif closure == 'py':
            U_mat = kBT * np.log(1 - c_mat/g_mat)
    if verbose:
        np.savez_compressed('calc-U-matrix.npz', r=r, k=k,
                            g_mat=g_mat, h_hat_mat=h_hat_mat,
                            G_minus_g_hat_mat=G_minus_g_hat_mat, rhos=rhos,
                            kBT=kBT, closure=closure, U_mat=U_mat)
    return U_mat


def calc_dU_newton(r, g_tgt, g_cur, G_minus_g, n, kBT, rho, cut_off,
                   closure, newton_mod, cut_jacobian, verbose=False):
    """
    Calculate a potential update dU using Newtons method.

    Supports symmetric molecules with n equal beads.

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
        cut_jacobian: Wether to cut the Jacobian. If False, then the full Jacobian will
            be used.

    Returns:
        The calculated potential update.

    """
    r0_removed, (r, g_tgt, g_cur, G_minus_g) = r0_removal(r, g_tgt, g_cur, G_minus_g)
    # calc slices and Delta_r
    nocore, crucial = calc_slices(r, g_tgt, g_cur, cut_off, verbose=verbose)
    # difference of rdf to target
    Delta_g = g_cur - g_tgt
    # FT of total correlation function 'h'
    _, h_hat = fourier(r, g_cur - 1)
    F = gen_fourier_matrix(r, fourier)
    # dc/dg
    if n == 1:
        # single bead case
        dcdg = np.linalg.inv(F) @ np.diag(1 / (1 + rho * h_hat)**2) @ F
    else:
        _, G_minus_g_hat = fourier(r, G_minus_g)
        dcdg = np.linalg.inv(F) @ np.diag(1 / (1 + n * rho * G_minus_g_hat
                                               + n * rho * h_hat)**2) @ F
    # calculate jacobian^-1
    # in the core where RDF=0, the jacobin will have -np.inf on the diagonal
    # numpy correctly inverts this to zero
    if closure == 'hnc':
        if newton_mod:
            with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
                jac_inv1 = kBT * (1 + np.log(g_tgt / g_cur) / Delta_g)
            jac_inv2 = -kBT * dcdg
            # Some fixes, because we want to define a jacobian matrix
            # Unmodified Newton is less awkward
            # Ensure this is zero, not nan, on the diagonal where Delta_g is zero
            jac_inv1[Delta_g == 0] = 0
            # Ensure this is -np.inf, not nan, on the diagonal in the core region
            jac_inv1[:nocore.start] = -np.inf
            jac_inv = np.diag(jac_inv1) + jac_inv2
        else:
            with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
                jac_inv = kBT * (np.diag(1 - 1 / g_cur) - dcdg)
    elif closure == 'py':
        raise NotImplementedError
    jac = np.linalg.inv(jac_inv)
    if cut_jacobian:
        # cut jacobian and transform back
        jac_cut = jac[crucial, crucial]
        jac_inv_cut = np.linalg.inv(jac_cut)
        dU = - (jac_inv_cut @ Delta_g[crucial])
    else:
        with np.errstate(invalid='ignore'):
            dU = - (jac_inv @ Delta_g)[crucial]
    if verbose:
        np.savez_compressed('newton-arrays.npz', jac=jac, jac_inv=jac_inv, dU=dU)
    # fill core and behind cut-off
    dU = np.concatenate((np.full(nocore.start, np.nan), dU,
                         np.full(len(r) - crucial.stop, np.nan)))
    if r0_removed:
        dU = np.concatenate(([np.nan], dU))
    return dU


def calc_dU_gauss_newton(r, g_tgt, g_cur, G_minus_g, n, kBT, rho,
                         cut_off, constraints,
                         verbose=False):
    """
    Calculate a potential update dU using the Gauss-Newton method.

    Constraints can be added.

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
    r0_removed, (r, g_tgt, g_cur, G_minus_g) = r0_removal(r, g_tgt, g_cur, G_minus_g)
    # calc slices and Delta_r
    nocore, crucial = calc_slices(r, g_tgt, g_cur, cut_off, verbose=verbose)
    Delta_r = calc_grid_spacing(r)
    # pair correlation function 'h'
    h = g_cur - 1
    # special Fourier of h
    _, h_hat = fourier(r, h)
    # Fourier matrix
    F = gen_fourier_matrix(r, fourier)
    # dc/dg
    if n == 1:
        # single bead case
        dcdg = np.linalg.inv(F) @ np.diag(1 / (1 + rho * h_hat)**2) @ F
    else:
        _, G_minus_g_hat = fourier(r, G_minus_g)
        dcdg = np.linalg.inv(F) @ np.diag(1 / (1 + n * rho * G_minus_g_hat
                                               + n * rho * h_hat)**2) @ F
    # jacobian^-1 (matrix U in Delbary et al., with respect to potential)
    with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
        jac_inv = kBT * (np.diag(1 - 1 / g_cur[nocore]) - dcdg[nocore, nocore])
    # A0 matrix
    A0 = Delta_r * np.triu(np.ones((len(r[nocore]), len(r[crucial])-1)), k=0)
    # Jacobian with respect to force
    J = np.linalg.inv(jac_inv) @ A0
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
    dU = A0 @ w
    # fill core with nans
    dU = np.concatenate((np.full(nocore.start, np.nan), dU))
    # dump files
    if verbose:
        np.savez_compressed('gauss-newton-arrays.npz', A=A, b=b, C=C, d=d,
                            jac_inv=jac_inv, A0=A0, J=J)
    if r0_removed:
        dU = np.concatenate(([np.nan], dU))
    return dU


def shift_U_cutoff_zero(dU, r, U, cut_off):
    """Make potential zero at and beyond cut-off."""
    dU_shift = dU.copy()
    # shift dU to be zero at cut_off and beyond
    ndx_co = find_after_cut_off_ndx(r, cut_off)
    U_before_cut_off = U[ndx_co-1] + dU[ndx_co-1]
    dU_shift -= U_before_cut_off
    dU_shift[ndx_co:] = -1 * U[ndx_co:]
    return dU_shift


def get_args(iie_args=None):
    """Define and parse command line arguments."""
    description = """
    Calculate U or ΔU with Integral Equations.
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
        pars.add_argument('--volume', type=float,
                          required=True,
                          metavar='VOL',
                          help='the volume of the box')
        pars.add_argument('--topol', type=argparse.FileType('r'),
                          required=True,
                          metavar='TOPOL',
                          help='XML topology file')
        pars.add_argument('--options', type=argparse.FileType('r'),
                          required=True,
                          metavar='SETTINGS',
                          help='XML settings file')
        pars.add_argument('--g-tgt-ext', type=str,
                          required=True,
                          metavar='RDF_TGT_EXT',
                          help='extension of RDF target files')
        pars.add_argument('--out-ext', type=str,
                          required=True,
                          metavar='U_OUT_EXT',
                          help='extension of U or ΔU output files')
        pars.add_argument('--g-extrap-factor', type=float, required=False,
                          help='factor by which to extrapolate RDFs')
        pars.add_argument('--g-min', type=float, required=True,
                          help='minimum value of RDF to be considered, determines valid'
                          'range')
    # intial potential guess subparsers
    for pars in [parser_pot_guess]:
        pars.add_argument('--g-tgt-intra-ext', type=str,
                          required=True,
                          metavar='RDF_TGT_INTRA_EXT',
                          help='extension of intramol. RDF target files')
    # update potential subparsers
    for pars in [parser_newton, parser_newton_mod, parser_gauss_newton]:
        pars.add_argument('--g-cur-ext', type=str,
                          required=True,
                          metavar='RDF_CUR_EXT',
                          help='extension of current RDF files')
        pars.add_argument('--g-cur-intra-ext', type=str,
                          required=True,
                          metavar='RDF_CUR_INTRA_EXT',
                          help='extension of current intramol. RDF files')
    # Newton's method only options
    for pars in [parser_newton, parser_newton_mod]:
        pars.add_argument('--cut-jacobian', dest='cut_jacobian', action='store_true',
                          help=('Cut the top-left part of the Jacobian before'
                                + ' multiplying with Δg.'))
        pars.set_defaults(cut_jacobian=False)
    # HNCGN only options
    parser_gauss_newton.add_argument('--pressure-constraint',
                                     dest='pressure_constraint',
                                     type=str, default=None)
    # parse
    if iie_args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(iie_args)
    # check for subcommand
    if args.subcommand is None:
        parser.print_help()
        raise Exception("subcommand needed")
    return args


def process_input(args):
    """Process arguments and perform some checks."""
    # TODO: g_cur, G_cur
    # args.options.read() can be called only once
    options = ET.fromstring(args.options.read())
    topology = ET.fromstring(args.topol.read())
    # get density_dict and n_intra_dict
    density_dict = get_density_dict(topology, args.volume)
    n_intra_dict = get_n_intra_dict(topology)
    # get non_bonded_dict
    non_bonded_dict = {nb_name: nb_ts for nb_name, nb_ts in get_non_bonded(options)}
    non_bonded_dict_inv = {v: k for k, v in non_bonded_dict.items()}
    if len(non_bonded_dict) != len(non_bonded_dict_inv):
        raise Exception("Some non-bonded name was not unique or some non-bonded "
                        "interactions had the same bead types.")
    # dict of table extensions
    table_infos = {
        'g_tgt': {'extension': args.g_tgt_ext, 'check-grid': True},
        'G_minus_g_tgt': {'extension': args.g_tgt_intra_ext, 'check-grid': True},
    }
    # load input arrays
    input_arrays = {}  # will hold all input data
    for table_name, table_info in table_infos.items():
        input_arrays[table_name] = {}
        for non_bonded_name in non_bonded_dict.keys():
            x, y, flag = readin_table(non_bonded_name + table_info['extension'])
            input_arrays[table_name][non_bonded_name] = {'x': x, 'y': y, 'flag': flag}
    # check for same grid and define r
    r = None
    for table_name, table_info in table_infos.items():
        for non_bonded_name in non_bonded_dict.keys():
            x = input_arrays[table_name][non_bonded_name]['x']
            if table_info['check-grid']:
                if r is None:
                    # set first r
                    r = x
                else:
                    # compare with first r
                    if not np.allclose(x, r):
                        raise RuntimeError("Grids of tables do not match")
    # check if starts at r = 0.0, if so: remove
    all_first_x = np.array([input_arrays[table_name][non_bonded_name]['x'][0]
                            for table_name, table_info in table_infos.items()
                            for non_bonded_name in non_bonded_dict.keys()])
    # if all r[0] = 0
    if np.allclose(all_first_x, np.zeros_like(all_first_x)):
        for table_name, table_info in table_infos.items():
            for non_bonded_name in non_bonded_dict.keys():
                for key in ('x', 'y', 'flag'):
                    input_arrays[table_name][non_bonded_name][key] = (
                        input_arrays[table_name][non_bonded_name][key][1:])
        r = r[1:]
        r0_removed = True
    # if they do not start with 0 but are all the same value
    elif np.allclose(all_first_x, all_first_x[0]):
        r0_removed = False
    else:
        raise Exception('either all or no input tables should start at r=0')
    # process input further
    rhos = gen_beadtype_property_array(density_dict, non_bonded_dict)
    n_intra = gen_beadtype_property_array(n_intra_dict, non_bonded_dict)
    # settings
    # copy some directly from args
    settings_to_copy = ('closure', 'verbose', 'out_ext', 'g_min')
    settings = {key: vars(args)[key] for key in settings_to_copy}
    settings['non-bonded-dict'] = non_bonded_dict
    # settings['densities'] = densities  # not sure if needed later
    settings['rhos'] = rhos
    settings['n_intra'] = n_intra
    settings['r0-removed'] = r0_removed
    # others from options xml
    settings['kBT'] = float(options.find("./inverse/kBT").text)
    if args.subcommand == 'potential_guess':
        settings['cut_off'] = float(
            options.find("./inverse/initial_guess/ie/cut_off").text)
    else:
        settings['cut_off'] = float(
            options.find("./inverse/iie/cut_off").text)
    return r, input_arrays, settings


def potential_guess(r, input_arrays, settings):
    """Do the potential guess."""
    # prepare matrices
    g_mat = gen_interaction_matrix(r, input_arrays['g_tgt'],
                                   settings['non-bonded-dict'])
    h_mat = g_mat - 1
    k, h_hat_mat = fourier_all(r, h_mat)
    G_minus_g_mat = gen_interaction_matrix(
        r, input_arrays['G_minus_g_tgt'], settings['non-bonded-dict'])
    _, G_minus_g_hat_mat = fourier_all(r, G_minus_g_mat)
    # perform actual math
    U1_mat = calc_U_matrix(r, k, g_mat, h_hat_mat, G_minus_g_hat_mat,
                           settings['rhos'], settings['n_intra'],
                           settings['kBT'], settings['closure'],
                           verbose=settings['verbose'])
    if settings['verbose']:
        np.savez_compressed('potential-guess-arrays.npz', r=r, k=k,
                            g_mat=g_mat, h_hat_mat=h_hat_mat,
                            G_minus_g_hat_mat=G_minus_g_hat_mat, rhos=settings['rhos'],
                            n_intra=settings['n_intra'], U1_mat=U1_mat)
    # extrapolate and save potentials
    for non_bonded_name, U_dict in gen_interaction_dict(
            r, U1_mat, settings['non-bonded-dict']).items():
        U = U_dict['y']
        U_flag = gen_flag_isfinite(U)
        if settings['r0-removed']:
            r_out = np.concatenate(([0.0], r))
            U = np.concatenate(([np.nan], U))
            U_flag = np.concatenate((['o'], U_flag))
        else:
            r_out = r
        # change NaN in the core region to first valid value
        U = extrapolate_dU_left_constant(U, U_flag)
        comment = "created by: {}".format(" ".join(sys.argv))
        fn = non_bonded_name + settings['out_ext']
        saveto_table(fn, r_out, U, U_flag, comment)


def newton_update(r, input_arrays, args):
    """Do the Newton update."""
    # further tests on input
    # if g_extrap_factor != 1.0 then it should hold that cut-off == r[-1]
    # this is basically HNCN (ex) vs HNCN (ed)
    if not np.isclose(args.g_extrap_factor, 1.0):
        if not np.isclose(args.cut_off, r[-1]):
            raise Exception('if g_extrap_factor is not equal 1.0, the cut-off '
                            'needs to be the same as the range of all RDFs for '
                            'Newton method.')
        # extrapolate RDFs
        if args.g_extrap_factor < 1:
            raise Exception('g_extrap_factor needs to be larger than 1.0')
        Delta_r = calc_grid_spacing(r)
        r_short = np.copy(r)
        r = np.arange(r[0], r[-1] * args.g_extrap_factor, Delta_r)
        g_tgt = extrapolate_g(r_short, r, input_arrays['g_tgt'][0]['y'],
                              input_arrays['G_minus_g'][0]['y'],
                              args.n_intra[0], args.densities[0],
                              verbose=args.verbose)
        g_cur = extrapolate_g(r_short, r, input_arrays['g_cur'][0]['y'],
                              input_arrays['G_minus_g'][0]['y'],
                              args.n_intra[0], args.densities[0],
                              verbose=args.verbose)
        if args.verbose:
            np.savez_compressed('extrapolated.npz', r=r, g_tgt=g_tgt, g_cur=g_cur)
        G_minus_g = np.concatenate((input_arrays['G_minus_g'][0]['y'],
                                    np.zeros(len(r)-len(r_short))))
        cut_off = r_short[-1]
    else:
        g_tgt = input_arrays['g_tgt'][0]['y']
        g_cur = input_arrays['g_cur'][0]['y']
        G_minus_g = input_arrays['G_minus_g'][0]['y']
        cut_off = args.cut_off
    dU1 = calc_dU_newton(r, g_tgt, g_cur, G_minus_g,
                         args.n_intra[0], args.kBT,
                         args.densities[0], cut_off, args.closure,
                         args.subcommand == 'newton-mod',
                         args.cut_jacobian,
                         verbose=args.verbose)
    # go back to short r range if we extrapolated
    if not np.isclose(args.g_extrap_factor, 1.0):
        dU1 = dU1[:len(r_short)]
        r = r_short
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
                                       cut_off)
    dU3 = U_temp - input_arrays['U_cur'][0]['y']
    comment = "created by: {}".format(" ".join(sys.argv))
    saveto_table(args.U_out[0], r, dU3, dU_flag4, comment)


def gauss_newton_update(r, input_arrays, args):
    """Do the Gauss-Newton update."""
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


def main():
    # get command line arguments
    args = get_args()

    # process and prepare input
    r, input_arrays, settings = process_input(args)

    # guess potential from distribution
    if args.subcommand == 'potential_guess':
        potential_guess(r, input_arrays, settings)

    # newton update
    if args.subcommand in ('newton', 'newton-mod'):
        newton_update(r, input_arrays, args)

    # gauss-newton update
    if args.subcommand == 'gauss-newton':
        gauss_newton_update(r, input_arrays, args)


if __name__ == '__main__':
    main()
