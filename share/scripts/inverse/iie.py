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


def fourier(r, f):
    """Compute the radially 3D FT and the frequency grid of a radially symmetric function.
    Some special extrapolations are used to make the results consistent. This function
    is isometric meaning it can be used to calculate the FT and the inverse FT.
    That means inputs can also be omega and f_hat which results in r and f.

    Args:
        r: Input grid. Must be evenly spaced. Can start at zero or at Δr, but nowhere
            else.
        f: Input function. Must have same length as r and correspond to its values.

    Returns:
        (omega, f_hat): The reciprocal grid and the FT of f.
    """
    Delta_r = calc_grid_spacing(r)
    r0_added = False
    if np.isclose(r[0], Delta_r):
        r = np.concatenate(([0], r))
        f = np.concatenate(([0], f))
        r0_added = True
    elif np.isclose(r[0], 0.0):
        pass
    else:
        raise Exception('this function can not handle this input')
    # if the input is even, np.fft.rfftfreq would end with the Nyquist frequency.
    # But there the imaginary part of the FT is always zero, so we alwas append a zero
    # to obtain a odd grid.
    if len(r) % 2 == 0:  # even
        r = np.concatenate((r, [r[-1]+Delta_r]))
        f = np.concatenate((f, [0]))
        n = (len(r)-1)*2-1
    else:  # odd
        n = len(r)*2-1
    omega = np.fft.rfftfreq(n=n, d=Delta_r)
    with np.errstate(divide='ignore', invalid='ignore'):
        f_hat = -2 / omega / 1 * Delta_r * np.imag(np.fft.rfft(r * f, n=n))
    if r0_added:
        f_hat = f_hat[1:]
        omega = omega[1:]
    return omega, f_hat


def gen_fourier_matrix(r, fourier_function):
    """make a fourier matrix"""
    fourier_matrix = np.identity(len(r))
    for col_index, col in enumerate(fourier_matrix.T):
        _, fourier_matrix.T[col_index] = fourier_function(r, col)
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


def r0_removal(*arrays):
    r0_removed = False
    if np.isclose(arrays[0][0], 0.0):
        r0_removed = True
        arrays = tuple(map(lambda a: a[1:], arrays))
    return r0_removed, arrays


def calc_c(r, g_tgt, G_minus_g, n, rho):
    r0_removed, (r, g_tgt, G_minus_g) = r0_removal(r, g_tgt, G_minus_g)
    # total correlation function h
    h = g_tgt - 1
    # FT of h
    omega, h_hat = fourier(r, h)
    # direct correlation function c from OZ
    if n == 1:
        # single bead case
        c_hat = h_hat / (1 + rho * h_hat)
    else:
        _, G_minus_g_hat = fourier(r, G_minus_g)
        c_hat = h_hat / ((1 + n * rho * G_minus_g_hat)**2
                         + (1 + n * rho * G_minus_g_hat) * n * rho * h_hat)
    _, c = fourier(omega, c_hat)
    if r0_removed:
        c = np.concatenate(([np.nan], c))
    return c


def calc_g(r, c, G_minus_g, n, rho):
    r0_removed, (r, c, G_minus_g) = r0_removal(r, c, G_minus_g)
    omega, c_hat = fourier(r, c)
    if n == 1:
        h_hat = c_hat / (1 - rho * c_hat)
    else:
        _, G_minus_g_hat = fourier(r, G_minus_g)
        h_hat = ((c_hat * (1 + n * rho * G_minus_g_hat)**2)
                 / (1 - n * rho * (1 + n * rho * G_minus_g_hat) * c_hat))
    _, h = fourier(omega, h_hat)
    g = h + 1
    if r0_removed:
        g = np.concatenate(([np.nan], g))
    return g


def calc_dc_ext(r_short, r_long, c_k_short, g_k_short, g_tgt_short, G_minus_g_short, n,
                rho):
    """Calculate Δc_ext with netwon method. Jacobian has an implicit extrapolation of c
    with zeros on twice its original range."""
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
    """Extrapolate an RDF with integral equation theory.
    Assumes c = 0 in the extrapolated region.

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
    else:
        return g_long_k[-1]


def calc_slices(r, g_tgt, g_cur, cut_off, verbose=False):
    # there are different regions in r used
    #              |       crucial     |                     # regions (slices)
    # |   core     |                 nocore               |
    # 0---------core_end------------cut_off-----------r[-1]  # distances
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
    h = g_tgt - 1
    c = calc_c(r, g_tgt, G_minus_g, n, rho)
    with np.errstate(divide='ignore', invalid='ignore'):
        if closure == 'hnc':
            U = kBT * (-np.log(g_tgt) + h - c)
        elif closure == 'py':
            U = kBT * np.log(1 - c/g_tgt)
    return U


def calc_dU_newton(r, g_tgt, g_cur, G_minus_g, n, kBT, rho, cut_off,
                   closure, newton_mod, cut_jacobian, verbose=False):
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
    omega, h_hat = fourier(r, g_cur - 1)
    F = gen_fourier_matrix(r, fourier)
    # dc/dg
    if n == 1:
        # single bead case
        dcdg = np.matmul(np.matmul(np.linalg.inv(F),
                                   np.diag(1 / (1 + rho * h_hat)**2)),
                         F)
    else:
        _, G_minus_g_hat = fourier(r, G_minus_g)
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
    jac = np.linalg.inv(jac_inv)
    if cut_jacobian:
        # cut jacobian and transform back
        jac_cut = jac[crucial, crucial]
        jac_inv_cut = np.linalg.inv(jac_cut)
        dU = - np.matmul(jac_inv_cut, Delta_g[crucial])
    else:
        with np.errstate(invalid='ignore'):
            dU = - np.matmul(jac_inv, Delta_g)[crucial]
    if verbose:
        np.savez_compressed('newton-arrays.npz', jac=jac, jac_inv=jac_inv, dU=dU)
    # fill core and behind cut-off
    dU = np.concatenate((np.full(nocore.start, np.nan), dU,
                         np.full(len(r) - crucial.stop, np.nan)))
    if r0_removed:
        dU = np.concatenate(([0], dU))
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
    r0_removed, (r, g_tgt, g_cur, G_minus_g) = r0_removal(r, g_tgt, g_cur, G_minus_g)
    # calc slices and Delta_r
    nocore, crucial = calc_slices(r, g_tgt, g_cur, cut_off, verbose=verbose)
    Delta_r = calc_grid_spacing(r)
    # pair correlation function 'h'
    h = g_cur - 1
    # special Fourier of h
    omega, h_hat = fourier(r, h)
    # Fourier matrix
    F = gen_fourier_matrix(r, fourier)
    # dc/dg
    if n == 1:
        # single bead case
        dcdg = np.matmul(np.matmul(np.linalg.inv(F),
                                   np.diag(1 / (1 + rho * h_hat)**2)),
                         F)
    else:
        _, G_minus_g_hat = fourier(r, G_minus_g)
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
    dU = np.concatenate((np.full(nocore.start, np.nan), dU))
    # dump files
    if verbose:
        np.savez_compressed('gauss-newton-arrays.npz', A=A, b=b, C=C, d=d,
                            jac_inv=jac_inv, A0=A0, J=J)
    if r0_removed:
        dU = np.concatenate(([0], dU))
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
        pars.add_argument('--g-extrap-factor', type=float, required=False,
                          help='factor by which to extrapolate RDFs')
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
