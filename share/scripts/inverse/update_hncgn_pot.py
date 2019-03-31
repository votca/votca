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

import argparse
import numpy as np
import sys


BAR_PER_MD_PRESSURE = 16.6053904


def readin_table(filename):
    table_dtype = {'names': ('x', 'y', 'y_flag'),
                   'formats': ('f', 'f', 'S1')}
    x, y, y_flag = np.loadtxt(filename, dtype=table_dtype, comments=['#', '@'], unpack=True)
    return x, y, y_flag


def saveto_table(filename, x, y, y_flag, comment=""):
    data = np.column_stack((x.T, y.T, y_flag.T))
    np.savetxt(filename, data, header=comment, fmt='%s')


def compare_grids(grid_a, grid_b):
    if np.any(grid_a - grid_b > 0.0001):
        print("Different grids!")
        sys.exit(1)


def fourier(r, y, omega):
    Delta_r = r[1] - r[0]
    y_hat = np.zeros_like(omega)
    np.seterr(divide='ignore', invalid='ignore', under='ignore')
    for i, omega_i in enumerate(omega):
        y_hat[i] = 2 / omega_i * Delta_r * np.sum(r * y * np.sin(2 * np.pi * omega_i * r))
    np.seterr(all='raise')
    return y_hat


def gen_fourier_matrix(r_nz, omega_nz, fourier_function):
    fourier_matrix = np.identity(len(omega_nz))
    for col_index, col in enumerate(fourier_matrix.T):
        fourier_matrix.T[col_index] = fourier_function(r_nz, col, omega_nz)
    return fourier_matrix


def gauss_newton_constrained(A, C, b, d):
    m,n = A.shape
    p,n = C.shape
    b.shape = (m)
    d.shape = (p)

    if p > 1:
        raise Exception("not implemented for p > 1")

    A_elim = A.copy()
    b_elim = b.copy()
    for i in range(p):

        pivot = np.argmax(abs(C[i]))  # find max value of C

        A_elim = A - np.ones_like(A) * A[:, pivot][:, np.newaxis] * C[i] / C[i, pivot]
        b_elim = b - A[:, pivot] * d[i] / C[i, pivot]
        A_elim = np.delete(A_elim, pivot, 1)

    if p == n:
        print("Warning: solution determined fully by constraints")
        x_elim = []
    else:
        #x_elim = np.linalg.solve(A_elim.T @ A_elim, A_elim.T @ b_elim)
        x_elim = np.linalg.solve(np.matmul(A_elim.T, A_elim), np.matmul(A_elim.T, b_elim))

    if p == 0:
        # no constraints
        x = x_elim
    else:
        #x_pivot = (d[i] - np.delete(C, pivot, 1) @ x_elim) / C[i, pivot]
        x_pivot = (d[i] - np.matmul(np.delete(C, pivot, 1), x_elim)) / C[i, pivot]
        x = np.insert(x_elim, pivot, x_pivot)

    return x


def calc_dpot_hncgn_core(r, rdf_current_g, rdf_target_g, kBT, density, g_min,
                         r_cut_off, pot_current_U, constraints, dump_steps=False):
    """Apply the (constrained) HNCGN method.
    Matrix multiplication was done with @ in earlier versions, now
    with numpy.matmul"""

    # for convenience grids (real and reciprocal) are without the zero value
    r = r[1:]
    delta_r = r[1] - r[0]
    g = rdf_current_g[1:]
    g_tgt = rdf_target_g[1:]
    u = pot_current_U[1:]

    # there are different regions in r used in the method
    #              |       main        |                     # regions
    # |   core     |                 nocore               |
    # 0 --------core_end------------cut_off----------len(r)  # those are indices
    #
    # Δ from the paper equals nocore
    # Δ' from the paper equals main
    # note: Vector w is in region main, but with one element less, because of
    #       the antiderivative operator A0

    core_end = np.where(g_tgt > g_min)[0][0]
    cut_off = np.searchsorted(r, r_cut_off)
    main = slice(core_end, cut_off+1)
    nocore = slice(core_end, len(r))
    #print("core_end:", core_end, r[core_end])
    #print("cut_off:", cut_off, r[cut_off])
    #print("r_end:", len(r)-1, r[-1])
    #print("r", 0, len(r), min(r), max(r))
    #print("main:", core_end, cut_off+1, min(r[main]), max(r[main]))
    #print("nocore:", core_end, len(r), min(r[nocore]), max(r[nocore]))

    # reciprocal grid up to Nyquist frequency
    delta_r = (r[-1] - r[0]) / (len(r) - 1)
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
    #T = np.linalg.inv(F) @ H @ F
    T = np.matmul(np.linalg.inv(F), np.matmul(H, F))

    # D matrix
    D = np.diag(g_tgt[core_end:])

    # U matrix
    U = -kBT * np.linalg.inv(D) + kBT * T[core_end:, core_end:]

    # A0 matrix
    A0 = delta_r * np.triu(np.ones((len(r[nocore]), len(r[main])-1)), k=0)

    # Jacobian matrix
    #J = np.linalg.inv(U) @ A0
    J = np.matmul(np.linalg.inv(U), A0)
    print('J', J.shape)

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
    #A = J.T @ J
    #b = J.T @ res[nocore]
    A = J
    b = res[nocore]
    print('d =', d)
    w = gauss_newton_constrained(A, C, b, d)

    # dU
    #dU = A0 @ w
    dU = np.matmul(A0, w)

    # fill core with nans
    dU = np.concatenate((np.full(core_end + 1, np.nan), dU))

    # dump files
    if dump_steps:
        np.savetxt("hncgn_h_hat.xvg", (omega, h_hat), header="omega, h_hat")

    return dU


def extrapolate_U_cubic(r, dpot_dU, dpot_flag, rdf_current_g, pot_current_U):
    """Extrapolate the potential in the core region up to
    the point, where the RDF becomes larger than 0.1 or where
    the new potential is convex. A cubic function is used.
    The first five valid points are uset for the fit.

    Fitting is done, because p-HNCGN often has artifacs at
    small RDF values, especially when pressure is far of

    Returns dU. U_{k+1} = U_k + dU is done by Votca."""
    # make copy
    dpot_dU_extrap = dpot_dU.copy()
    # region to fit
    fit_start1 = np.where(rdf_current_g > 1e-1)[0][0]
    fit_start2 = np.where(np.nan_to_num(np.diff(np.diff(
        pot_current_U + dpot_dU))) > 0)[0][0]
    fit_start = max(fit_start1, fit_start2)
    fit_end = fit_start + 5
    fit_region = slice(fit_start, fit_end)
    fit_x = r[fit_region]
    fit_y = (pot_current_U + dpot_dU)[fit_region]
    # fit p[0] x² + p[1] x + p[2]
    p = np.polyfit(fit_x, fit_y, 3)
    U_extrap = np.polyval(p, r)
    # region to extrapolate
    extrap_region = slice(0, fit_start)
    # extrapolate
    dpot_dU_extrap[extrap_region] = (U_extrap[extrap_region]
                                     - pot_current_U[extrap_region])
    return dpot_dU_extrap


def extrapolate_U_times_g(r, dpot_dU, dpot_flag, rdf_current_g):
    """Extrapolate the potential in the core region up to
    the point, where the RDF becomes larger than 0.1.

    Returns dU. U_{k+1} = U_k + dU is done by Votca."""
    # fill core with first value
    first_dU_index = np.asarray(dpot_flag == 'i').nonzero()[0][0]
    first_dU = dpot_dU[first_dU_index]
    # replace out of range dU values with constant first value
    dpot_dU_extrap = np.where(dpot_flag == 'i', dpot_dU, first_dU)
    # region to extrapolate
    extrap_end = np.where(rdf_current_g > 1e-0)[0][0]
    extrap_region = slice(0, extrap_end)
    # extrapolate
    dpot_dU_extrap[extrap_region] *= rdf_current_g[extrap_region]
    dpot_dU_extrap = np.nan_to_num(dpot_dU_extrap)
    return dpot_dU_extrap


def calc_dpot_hncgn(r, rdf_target_g, rdf_target_flag,
                    rdf_current_g, rdf_current_flag,
                    pot_current_U, pot_current_flag,
                    kBT, density, cut_off, g_min,
                    constraints, extrapolation):
    """calculate dU for the hncgn method with constraint pressure"""

    # allways raise an error
    np.seterr(all='raise')

    # prepare dpot
    dpot_dU = np.zeros_like(pot_current_U)
    dpot_flag = np.array([''] * len(dpot_dU))

    # full range dU
    dU_full = calc_dpot_hncgn_core(r, rdf_current_g, rdf_target_g, kBT,
                                   density, g_min, cut_off, pot_current_U,
                                   constraints)

    # calculate dpot
    for i in range(len(r)):
        if (np.isnan(dU_full[i])) or ('u' in str(pot_current_flag[i])):
            dpot_dU[i] = np.nan
            dpot_flag[i] = 'o'
        else:
            dpot_dU[i] = dU_full[i]
            dpot_flag[i] = 'i'

    # extrapolation in and near core region
    if extrapolation == 'none':
        dpot_dU_extrap = dpot_dU
    elif extrapolation == 'cubic':
        dpot_dU_extrap = extrapolate_U_cubic(r, dpot_dU, dpot_flag,
                                             rdf_current_g, pot_current_U)
    elif extrapolation == 'times_g':
        dpot_dU_extrap = extrapolate_U_times_g(r, dpot_dU, dpot_flag,
                                               rdf_current_g)
    else:
        raise Exception("unknow extrapolation scheme for inside and near core"
                        "region:" + extrapolation)
    return dpot_dU_extrap, dpot_flag


description = """\
This script calculatess dU with the HNCGN scheme.
It uses some magic tricks:
- beyond the cut_off dU is set to -U such that U becomes zero.
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('rdf_target', type=argparse.FileType('r'))
parser.add_argument('rdf_current', type=argparse.FileType('r'))
parser.add_argument('pot_current', type=argparse.FileType('r'))
parser.add_argument('dpot', type=argparse.FileType('wb'))
parser.add_argument('kBT', type=float)
parser.add_argument('density', type=float)
parser.add_argument('cut_off', type=float)
parser.add_argument('--pressure-constraint', dest='pressure_constraint',
                    type=str, default=None)
parser.add_argument('--extrap-near-core', dest='extrap_near_core',
                    type=str, default='times_g',
                    choices=['times_g', 'none', 'cubic'])

if __name__ == '__main__':
    args = parser.parse_args()

    # load rdf and potential
    rdf_target_r, rdf_target_g, rdf_target_flag = readin_table(args.rdf_target)
    rdf_current_r, rdf_current_g, rdf_current_flag = readin_table(args.rdf_current)
    pot_current_r, pot_current_U, pot_current_flag = readin_table(args.pot_current)

    # sanity checks on grid
    compare_grids(rdf_target_r, rdf_current_r)
    compare_grids(rdf_target_r, pot_current_r)
    r = rdf_target_r

    constraints = []
    if args.pressure_constraint is not None:
        p_target = float(args.pressure_constraint.split(',')[0])
        p_current = float(args.pressure_constraint.split(',')[1])
        constraints.append({'type': 'pressure', 'target': p_target,
                            'current': p_current})

    # calculate dpot
    dpot_dU, dpot_flag = calc_dpot_hncgn(r, rdf_target_g, rdf_target_flag,
                                         rdf_current_g, rdf_current_flag,
                                         pot_current_U, pot_current_flag,
                                         args.kBT, args.density, args.cut_off,
                                         1e-10, constraints,
                                         args.extrap_near_core)

    # save dpot
    comment = "created by: {}".format(" ".join(sys.argv))
    saveto_table(args.dpot, r, dpot_dU, dpot_flag, comment)
