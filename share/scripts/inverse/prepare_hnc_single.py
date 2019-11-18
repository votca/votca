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


def fourier(r, y, omega):
    Delta_r = r[1] - r[0]
    y_hat = np.zeros_like(omega)
    np.seterr(divide='ignore', invalid='ignore', under='ignore')
    for i, omega_i in enumerate(omega):
        y_hat[i] = 2 / omega_i * Delta_r * np.sum(r * y * np.sin(2 * np.pi * omega_i * r))
    np.seterr(all='raise')
    return y_hat


def calc_pot_hnc_core(r, rdf_current_g, kBT, density, dump_steps=False):
    """calculates U for all r (full width)"""

    # density 'ρ0'
    rho = density

    # pair correlation function 'h'
    h = rdf_current_g - 1

    # special Fourier of h
    omega = np.arange(1, len(r)) / (2 * max(r))
    h_hat = fourier(r, h, omega)

    # special Fourier of c (direct correlation function)
    c_hat = h_hat / (1 + rho * h_hat)

    # y
    y = h - fourier(omega, c_hat, r)

    # zero order U (similar IBI)
    np.seterr(divide='ignore', invalid='ignore')
    U_order_0 = -kBT * np.log(rdf_current_g)
    np.seterr(all='raise')

    # first order U
    U_order_1 = kBT * y

    # total U
    U = U_order_0 + U_order_1

    # dump files
    if dump_steps:
        np.savetxt("hnc_h_hat.xvg", (omega, h_hat), header="omega, h_hat")
        np.savetxt("hnc_c_hat.xvg", (omega, c_hat), header="omega, c_hat")

    return U


def calc_pot_hnc(r, rdf_target_g, rdf_target_flag,
                 kBT, density, cut_off, g_min):
    # allways raise an error
    np.seterr(all='raise')

    # prepare dpot
    pot_U = np.zeros_like(rdf_target_g)
    pot_flag = np.array([''] * len(rdf_target_g))

    # full range U
    U_full = calc_pot_hnc_core(r, rdf_target_g, kBT, density)

    # calculate dpot
    for i in range(len(r)):
        if rdf_target_g[i] > g_min:
            pot_U[i] = U_full[i]
            pot_flag[i] = 'i'
        else:
            pot_U[i] = np.nan
            pot_flag[i] = 'o'

    # find first valid U value
    first_U_index = np.where(pot_flag == 'i')[0][0]
    first_U = pot_U[first_U_index]

    # replace out of range U values
    pot_U = np.where(pot_flag == 'i', pot_U, first_U)

    # shift U to be zero at cut_off and beyond
    index_cut_off = np.searchsorted(r, cut_off)
    U_cut_off = pot_U[index_cut_off]
    pot_U -= U_cut_off
    pot_U[index_cut_off:] = 0

    return pot_U, pot_flag


description = """\
This script calculatess a guess for U with using the hypernetted chain (HNC)
approximation.
"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument('rdf_target', type=argparse.FileType('r'))
parser.add_argument('kBT', type=float)
parser.add_argument('density', type=float)
parser.add_argument('cut_off', type=float)
parser.add_argument('pot', type=argparse.FileType('wb'))  # output

if __name__ == '__main__':
    args = parser.parse_args()

    # load rdf and potential
    r, rdf_target_g, rdf_target_flag = readin_table(args.rdf_target)


    # calculate pot
    pot_U, pot_flag = calc_pot_hnc(r, rdf_target_g, rdf_target_flag,
                                   args.kBT, args.density, args.cut_off, 1e-10)

    # save pot
    comment = "created by: {}".format(" ".join(sys.argv))
    saveto_table(args.pot, r, pot_U, pot_flag, comment)
