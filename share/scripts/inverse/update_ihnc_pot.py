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
    np.savetxt(filename, (x, y, y_flag), header=comment, fmt='%s')


def compare_grids(grid_a, grid_b):
    if np.any(grid_a - grid_b > 0.0001):
        print("Different grids!")
        sys.exit(1)


def fft_paper(r, y):
    r = r[1:]
    y = y[1:]

    m = len(r)
    Delta_r = r[1] - r[0]
    l_range = range(-m, m+1, 1)
    j_range = l_range

    omega = np.array(l_range) / (m + 1) * 1 / (2 * Delta_r)

    r_pm = np.concatenate((-r[::-1], [0], r))
    y_pm = np.concatenate((y[::-1], [0], y))
    psi_pm = r_pm * y_pm

    y_hat = np.zeros_like(omega, dtype=np.complex128)

    for l_index, l in enumerate(l_range):
        omega_l = omega[l_index]
        for j_index, j in enumerate(j_range):
            y_hat[l_index] += psi_pm[j_index] * np.exp(-2 * np.pi * 1j * omega_l * r_pm[j_index])
        np.seterr(divide='ignore', invalid='ignore')
        y_hat[l_index] *= Delta_r / (1j * omega_l)
        np.seterr(all='raise')

    return omega, y_hat


def ifft_paper(omega, y_hat):

    zero_frequency_index = (len(omega) - 1) // 2
    y_hat[zero_frequency_index] = 0

    Delta_omega = omega[1] - omega[0]
    m = (len(omega) - 1) // 2
    l_range = range(-m, m+1, 1)
    j_range = l_range

    r_pm = np.array(l_range) / (m + 1) * 1 / (2 * Delta_omega)
    psi_pm = omega * y_hat

    y = np.zeros_like(r_pm, dtype=np.complex128)

    for l_index, l in enumerate(l_range):
        r_pm_l = r_pm[l_index]
        for j_index, j in enumerate(j_range):
            y[l_index] += psi_pm[j_index] * np.exp(-2 * np.pi * 1j * r_pm_l * omega[j_index])
        np.seterr(divide='ignore', invalid='ignore')
        y[l_index] *= Delta_omega / (1j * r_pm_l)
        np.seterr(all='raise')

    return r_pm[zero_frequency_index:], y[zero_frequency_index:]


def calc_dpot_ihnc_core(r, rdf_current_g, rdf_target_g, kBT, density):
    """calculates dU for all r (full width)"""

    # density 'ρ0'
    rho = density

    # difference of rdf to target 'f'
    f = rdf_target_g - rdf_current_g

    # pair correlation function 'h'
    h = rdf_current_g - 1

    # special Fourier of h
    omega, h_hat = fft_paper(r, h)

    # special Fourier of f
    omega, f_hat = fft_paper(r, f)

    # special Fourier of φ
    np.seterr(divide='ignore', invalid='ignore')
    phi_hat = (2 + rho * h_hat) / (1 + rho * h_hat)**2 * rho * h_hat * f_hat
    np.seterr(all='raise')

    # φ
    r_, phi = ifft_paper(omega, phi_hat)

    # zero order dU (similar IBI)
    np.seterr(divide='ignore', invalid='ignore')
    dU_order_0 = kBT * np.log(rdf_current_g / rdf_target_g)
    np.seterr(all='raise')

    # first order dU
    dU_order_1 = kBT * phi

    #print(dU_order_0, dU_order_1)

    dU = dU_order_0 + dU_order_1
    return dU


def calc_dpot_ihnc(r, rdf_target_g, rdf_target_flag,
                   rdf_current_g, rdf_current_flag,
                   pot_current_U, pot_current_flag,
                   kBT, density, g_min):
    # allways raise an error
    np.seterr(all='raise')

    # prepare dpot
    dpot_dU = np.zeros_like(pot_current_U)
    dpot_flag = np.array([''] * len(dpot_dU))

    # full range dU
    dU_full = calc_dpot_ihnc_core(r, rdf_current_g, rdf_target_g, kBT, density)

    # calculate dpot
    for i in range(len(r)):
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

    return dpot_dU, dpot_flag


description = """\
This script calculatess dU out of two rdfs with the rules of inverse boltzmann.
In addition, it does some magic tricks:
- do not update if one of the two rdf is undefined
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('rdf_target', type=argparse.FileType('r'))
parser.add_argument('rdf_current', type=argparse.FileType('r'))
parser.add_argument('pot_current', type=argparse.FileType('r'))
parser.add_argument('dpot', type=argparse.FileType('wb'))
parser.add_argument('kBT', type=float)

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
    dpot_dU, dpot_flag = calc_dpot_ihnc(r, rdf_target_g, rdf_target_flag,
                                        rdf_current_g, rdf_current_flag,
                                        pot_current_U, pot_current_flag,
                                        args.kBT, 33.3, 1e-10)

    # save dpot
    comment = "created by: {}".format(" ".join(sys.argv))
    saveto_table(args.dpot, r, dpot_dU, dpot_flag, comment)
