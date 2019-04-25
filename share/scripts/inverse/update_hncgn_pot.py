#!/usr/bin/env python3
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

# suffixes:
# _cur: current or of step k if currently doing iteration k
# _tgt: target
# _ce: core_end (where RDF gets > 0)
# _co: cut_off
#
# prefixes:
# ndx_: index

import argparse
import numpy as np
import sys


BAR_PER_MD_PRESSURE = 16.6053904


def readin_table(filename):
    """read in votca table"""
    table_dtype = {'names': ('x', 'y', 'y_flag'),
                   'formats': ('f', 'f', 'S1')}
    x, y, y_flag = np.loadtxt(filename, dtype=table_dtype, comments=['#', '@'],
                              unpack=True)
    return x, y, y_flag


def saveto_table(filename, x, y, y_flag, comment=""):
    """save votca table"""
    data = np.column_stack((x.T, y.T, y_flag.T))
    np.savetxt(filename, data, header=comment, fmt='%s')


def compare_grids(grid_a, grid_b):
    """check two grids for no point differing more than 0.1 picometer"""
    if np.any(grid_a - grid_b > 0.0001):
        raise Exception("Different grids!")


def fourier(r, y, omega):
    """fourier transform real data y on grid r into reciprocal space omega"""
    Delta_r = r[1] - r[0]
    y_hat = np.zeros_like(omega)
    np.seterr(divide='ignore', invalid='ignore', under='ignore')
    for i, omega_i in enumerate(omega):
        y_hat[i] = (2 / omega_i * Delta_r
                    * np.sum(r * y * np.sin(2 * np.pi * omega_i * r)))
    np.seterr(all='raise')
    return y_hat


def gen_fourier_matrix(r_nz, omega_nz, fourier_function):
    """make a fourier matrix"""
    fourier_matrix = np.identity(len(omega_nz))
    for col_index, col in enumerate(fourier_matrix.T):
        fourier_matrix.T[col_index] = fourier_function(r_nz, col, omega_nz)
    return fourier_matrix


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
        print("Warning: solution determined fully by constraints")
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


def find_nearest_ndx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def calc_dU_hncgn(r, g_cur, g_tgt, kBT, density,
                  g_min, cut_off, constraints,
                  verbose):
    """Apply the (constrained) HNCGN method."""

    # for convenience grids (real and reciprocal) are without the zero value
    r = r[1:]
    g = g_cur[1:]
    g_tgt = g_tgt[1:]
    # there are different regions in r used in the method
    #              |       main        |                     # regions
    # |   core     |                 nocore               |
    #  ---------core_end------------cut_off-----------r[-1]  # distances
    # 0----------ndx_ce--------------ndx_co----------len(r)  # indices
    #
    # Δ from the paper equals nocore
    # Δ' from the paper equals main
    # note: Vector w is in region main, but with one element less, because of
    #       the antiderivative operator A0
    ndx_ce = np.where(g_tgt > g_min)[0][0]
    ndx_co = find_nearest_ndx(r, cut_off)
    main = slice(ndx_ce, ndx_co+1)
    nocore = slice(ndx_ce, len(r))
    delta_r = (r[-1] - r[0]) / (len(r) - 1)
    if verbose:
        print("core_end:", ndx_ce, r[ndx_ce])
        print("cut_off:", ndx_co, r[ndx_co])
        print("r_end:", len(r)-1, r[-1])
        print("r", 0, len(r), min(r), max(r))
        print("main:", ndx_ce, ndx_co+1, min(r[main]), max(r[main]))
        print("nocore:", ndx_ce, len(r), min(r[nocore]), max(r[nocore]))

    # reciprocal grid up to Nyquist frequency
    nyquist_freq = 1 / delta_r / 2
    omega = np.linspace(nyquist_freq/len(r), nyquist_freq, len(r))
    # This ensures max(omega) becomes not to large due to numerical accuracy
    # because if it does there are artifacts
    omega *= 0.9999
    # Fourier matrix
    F = gen_fourier_matrix(r, omega, fourier)
    # density 'ρ0'
    rho = density
    # pair correlation function 'h'
    h = g - 1
    # special Fourier of h
    h_hat = fourier(r, h, omega)
    # H matrix
    H = np.diag((2 + rho * h_hat) / (1 + rho * h_hat)**2 * rho * h_hat)
    # T matrix
    T = np.matmul(np.linalg.inv(F), np.matmul(H, F))
    # D matrix
    D = np.diag(g_tgt[ndx_ce:])
    # U matrix
    U = -kBT * np.linalg.inv(D) + kBT * T[ndx_ce:, ndx_ce:]
    # A0 matrix
    A0 = delta_r * np.triu(np.ones((len(r[nocore]), len(r[main])-1)), k=0)
    # Jacobian matrix
    J = np.matmul(np.linalg.inv(U), A0)
    # constraint matrix and vector
    C = np.zeros((len(constraints), len(r[main])-1))
    d = np.zeros(len(constraints))
    # build constraint matrix and vector from constraints
    for c, constraint in enumerate(constraints):
        if constraint['type'] == 'pressure':
            # current pressure
            p = constraint['current'] / BAR_PER_MD_PRESSURE
            # target pressure
            p_tgt = constraint['target'] / BAR_PER_MD_PRESSURE
            # g_tgt(r_{i+1})
            g_tgt_ip1 = g_tgt[main][1:]
            # g_tgt(r_{i})
            g_tgt_i = g_tgt[main][:-1]
            # r_{i+1}
            r_ip1 = r[main][1:]
            # r_{i}
            r_i = r[main][:-1]
            # l vector
            l = (g_tgt_i + g_tgt_ip1) * (r_ip1**4 - r_i**4)
            l *= 1/12 * np.pi * rho**2
            # set C row and d element
            C[c, :] = l
            d[c] = p_tgt - p
        else:
            raise Exception("not implemented constraint type")
    # residuum vector
    res = g_tgt - g
    # switching to notation of Gander et al. for solving
    A = J
    b = res[nocore]
    w = gauss_newton_constrained(A, C, b, d)
    # dU
    dU = np.matmul(A0, w)
    # fill core with nans
    dU = np.concatenate((np.full(ndx_ce + 1, np.nan), dU))
    # dump files
    if verbose:
        np.savez_compressed('hncgn-arrays.npz', A=A, b=b, C=C, d=d,
                            H=H, T=T, D=D, U=U, A0=A0, J=J)
    return dU


def extrapolate_U_constant(dU, dU_flag):
    dU_extrap = dU.copy()
    # find first valid dU value
    first_dU_index = np.where(dU_flag == 'i')[0][0]
    first_dU = dU[first_dU_index]

    # replace out of range dU values with constant first value
    dU_extrap = np.where(dU_flag == 'i', dU, first_dU)
    return dU_extrap


def extrapolate_U_power(r, dU, g, U, g_tgt, g_min, kBT, verbose):
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
    # r without zero at start
    with np.errstate(divide='ignore', over='ignore'):
        pmf_fit = np.nan_to_num(a * r**b - pmf_shift)

    # region to extrapolate
    ndx_ex1 = ndx_ce + 1
    ndx_ex2 = np.where(np.nan_to_num(np.diff(np.diff(U + dU))) > 0)[0][0]
    ndx_ex = max(ndx_ex1, ndx_ex2)
    # extrapolate
    U_extrap = U + dU
    U_extrap[:ndx_ex] = pmf_fit[:ndx_ex] + (U_extrap[ndx_ex] - pmf_fit[ndx_ex])
    dU_extrap = U_extrap - U
    return dU_extrap


def shift_U_cutoff_zero(dU, r, U, cut_off):
    """Make potential zero at and beyond cut-off"""
    dU_shift = dU.copy()
    # shift dU to be zero at cut_off and beyond
    ndx_co = find_nearest_ndx(r, cut_off)
    U_cut_off = U[ndx_co] + dU[ndx_co]
    dU_shift -= U_cut_off
    dU_shift[ndx_co:] = -1 * U[ndx_co:]
    return dU_shift


def fix_U_near_cut_off_full(dU, r, U, cut_off):
    """Modify the potential close to the cut-off in
    a way, such that it is more smooth. The derivative
    of the potential between the last two points will
    be equal to the derivative between the two points
    before. The original last two points of dU are
    therefore ignored.

    This also helps agains an artifact of p-HNCGN,
    where the last value of dU is a spike."""
    dU_fixed = dU.copy()
    ndx_co = find_nearest_ndx(r, cut_off)
    U_new = U + dU
    second_last_deriv = U_new[ndx_co-1] - U_new[ndx_co-2]
    shift = -1.0 * second_last_deriv - U_new[ndx_co-1]
    # modify up to second last value
    dU_fixed[:ndx_co] += shift
    return dU_fixed


description = """\
This script calculatess dU with the HNCGN scheme.
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('g_tgt', type=argparse.FileType('r'))
parser.add_argument('g_cur', type=argparse.FileType('r'))
parser.add_argument('U_cur', type=argparse.FileType('r'))
parser.add_argument('dU', type=argparse.FileType('wb'))
parser.add_argument('kBT', type=float)
parser.add_argument('density', type=float)
parser.add_argument('cut_off', type=float)
parser.add_argument('--pressure-constraint', dest='pressure_constraint',
                    type=str, default=None)
parser.add_argument('--extrap-near-core', dest='extrap_near_core',
                    type=str, choices=['none', 'power'])
parser.add_argument('--fix-near-cut-off', dest='fix_near_cut_off',
                    type=str, choices=['none', 'full-deriv'])
parser.add_argument('-v', '--verbose', dest='verbose',
                    action='store_const', const=True, default=False)

if __name__ == '__main__':
    args = parser.parse_args()

    # load rdf and potential
    g_tgt_r, g_tgt, g_tgt_flag = readin_table(args.g_tgt)
    g_cur_r, g_cur, g_cur_flag = readin_table(args.g_cur)
    U_cur_r, U_cur, U_cur_flag = readin_table(args.U_cur)

    np.seterr(all='raise')
    G_MIN = 1e-10

    # sanity checks on grid
    compare_grids(g_tgt_r, g_cur_r)
    compare_grids(g_tgt_r, U_cur_r)
    r = g_tgt_r

    # parse constraints
    constraints = []
    if args.pressure_constraint is not None:
        p_target = float(args.pressure_constraint.split(',')[0])
        p_current = float(args.pressure_constraint.split(',')[1])
        constraints.append({'type': 'pressure', 'target': p_target,
                            'current': p_current})

    # calc dU_pure
    dU_pure = calc_dU_hncgn(r, g_cur, g_tgt, args.kBT,
                            args.density, G_MIN, args.cut_off,
                            constraints, args.verbose)

    # set dU_flag to 'o' inside the core
    dU_flag = np.where(np.isnan(dU_pure), 'o', 'i')

    # select extrapolation
    if args.extrap_near_core == 'none':
        dU_extrap = np.nan_to_num(dU_pure)
    if args.extrap_near_core == 'constant':
        dU_extrap = extrapolate_U_constant(dU_pure, dU_flag)
    elif args.extrap_near_core == 'power':
        dU_extrap = extrapolate_U_power(r, dU_pure, g_cur, U_cur, g_tgt, G_MIN,
                                        args.kBT, args.verbose)
    else:
        raise Exception("unknow extrapolation scheme for inside and near core"
                        "region:" + args.extrap_near_core)
    # shifts to correct potential after cut-off
    dU_shift = shift_U_cutoff_zero(dU_extrap, r, U_cur, args.cut_off)
    # shifts to correct potential near cut-off
    if args.fix_near_cut_off == 'none':
        dU = dU_shift.copy()
    elif args.fix_near_cut_off == 'full-deriv':
        dU = fix_U_near_cut_off_full(dU_shift, r, U_cur, args.cut_off)
    else:
        raise Exception("unknow fix scheme for near cut-off: "
                        + args.fix_near_cut_off)

    if args.verbose:
        np.savez_compressed('hncgn-dU.npz', r=r, dU_pure=dU_pure,
                            dU_extrap=dU_extrap, dU_shift=dU_shift)

    # save dU
    comment = "created by: {}".format(" ".join(sys.argv))
    saveto_table(args.dU, r, dU, dU_flag, comment)
