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


def calc_dpot_hncgn_full(r, rdf_current_g, rdf_target_g, kBT, density, g_min, cut_off, dump_steps=False):
    """calculates dU for all r with the Gauss-Newton approach"""

    # for convenience grids (real and reciprocal) are without the zero value
    r = r[1:]
    delta_r = r[1] - r[0]
    g = rdf_current_g[1:]
    g_tgt = rdf_target_g[1:]

    # reciprocal grid up to Nyquist frequency
    delta_r = (r[-1] - r[0]) / (len(r) - 1)
    nyquist_freq = 1 / delta_r / 2
    omega = np.linspace(nyquist_freq/len(r), nyquist_freq, len(r))
    # This ensures max(omega) becomes not to large due to numerical accuracy
    # because if it does there are artifacts
    omega *= 0.99999

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
    T = np.linalg.inv(F) @ H @ F

    # real grid without core region
    core_end = np.where(g_tgt > g_min)[0][0]
    r_nocore = r[core_end:]

    # D matrix
    D = np.diag(g_tgt[core_end:])

    # U matrix
    U = -kBT * np.linalg.inv(D) + kBT * T[core_end:, core_end:]

    # A0 matrix
    index_cut_off = np.searchsorted(r_nocore, cut_off)
    r_nocore_cut = r_nocore[:index_cut_off]
    A0 = delta_r * np.triu(np.ones((len(r_nocore), len(r_nocore_cut))), k=0)

    # Jacobian matrix (really that simple?)
    J = - np.linalg.inv(U) @ A0

    # residuum vector
    res = g_tgt - g
    res_nocore = res[core_end:]

    # w
    # (J.T @ J) w == - J.T @ res_nocore
    #     a     x ==       b
    w = np.linalg.solve(J.T @ J, -J.T @ res_nocore)

    # dU
    print(w)
    dU = A0 @ w

    # fill core with nans
    dU = np.concatenate((np.full(core_end + 1, np.nan), dU))

    # dump files
    if dump_steps:
        pass

    return dU


def calc_dpot_hncgn(r, rdf_target_g, rdf_target_flag,
                      rdf_current_g, rdf_current_flag,
                      pot_current_U, pot_current_flag,
                      kBT, density, cut_off, g_min):
    # allways raise an error
    np.seterr(all='raise')

    # prepare dpot
    dpot_dU = np.zeros_like(pot_current_U)
    dpot_flag = np.array([''] * len(dpot_dU))

    # full range dU
    dU_full = calc_dpot_hncgn_full(r, rdf_current_g, rdf_target_g, kBT, density, g_min, cut_off)

    # calculate dpot
    for i in range(len(r)):

        # test
        #dpot_dU[i] = dU_full[i]
        #continue

        if rdf_target_g[i] > g_min and rdf_current_g[i] > g_min:
            dpot_dU[i] = dU_full[i]
            dpot_flag[i] = 'i'
        else:
            dpot_dU[i] = np.nan
            dpot_flag[i] = 'o'
        # check for unset value in current potential
        if 'u' in str(pot_current_flag[i]):
            dpot_dU[i] = np.nan
            dpot_flag[i] = 'o'

    # find first valid dU value
    first_dU_index = np.where(dpot_flag == 'i')[0][0]
    first_dU = dpot_dU[first_dU_index]

    # replace out of range dU values
    dpot_dU = np.where(dpot_flag == 'i', dpot_dU, first_dU)

    # shift dU to be zero at cut_off and beyond
    index_cut_off = np.searchsorted(r, cut_off)
    U_cut_off = pot_current_U[index_cut_off] + dpot_dU[index_cut_off]
    dpot_dU -= U_cut_off
    dpot_dU[index_cut_off:] = - pot_current_U[index_cut_off:]

    return dpot_dU, dpot_flag


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

if __name__ == '__main__':
    args = parser.parse_args()

    # load rdf and potential
    rdf_target_r, rdf_target_g, rdf_target_flag = readin_table(args.rdf_target)
    rdf_current_r, rdf_current_g, rdf_current_flag = readin_table(args.rdf_current)
    pot_current_r, pot_current_U, pot_current_flag = readin_table(args.pot_current)

    # sanity checks on grid
    compare_grids(rdf_target_r, rdf_current_r)
    r = rdf_target_r

    # calculate dpot
    dpot_dU, dpot_flag = calc_dpot_hncgn(r, rdf_target_g, rdf_target_flag,
                                         rdf_current_g, rdf_current_flag,
                                         pot_current_U, pot_current_flag,
                                         args.kBT, args.density, args.cut_off, 1e-10)

    # save dpot
    comment = "created by: {}".format(" ".join(sys.argv))
    saveto_table(args.dpot, r, dpot_dU, dpot_flag, comment)
