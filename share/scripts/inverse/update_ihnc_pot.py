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
    data = np.loadtxt(filename, dtype=np.str, comments=['#', '@'])
    x = data[:, 0].astype(np.float)
    y = data[:, 1].astype(np.float)
    y_flag = data[:, 2].astype('S1')
    return x, y, y_flag


def saveto_table(filename, x, y, y_flag, comment=""):
    data = np.column_stack((x.T, y.T, y_flag.T))
    format_string = '%s'
    np.savetxt(filename, data, header=comment, fmt=format_string)


def compare_grids(grid_a, grid_b):
    if np.any(grid_a - grid_b > 0.0001):
        print("Different grids!")
        sys.exit(1)


def calculate_dpot(r, rdf_target_g, rdf_target_flag,
                   rdf_current_g, rdf_current_flag,
                   pot_current_U, pot_current_flag,
                   kBT):
    # prepare dpot
    dpot_dU = np.zeros_like(pot_current_U)
    dpot_flag = np.array([''] * len(dpot_dU))

    # calculate dpot
    for i in range(len(r)):
        if rdf_target_g[i] > 1e-10 and rdf_current_g[i] > 1e-10:
            dpot_dU[i] = np.log(rdf_current_g[i] / rdf_target_g[i]) * kBT
            dpot_flag[i] = 'i'
        else:
            dpot_dU[i] = np.nan
            dpot_flag[i] = 'o'
        # check for unset value in current potential
        if 'u' in str(pot_current_flag[i]):
            dpot_dU[i] = np.nan
            dpot_flag[i] = 'o'

    # find first in range dU value
    first_dU_index = np.where(dpot_flag == 'i')[0][0]
    first_dU = dpot_dU[first_dU_index]
    first_dU_array = np.ones_like(dpot_dU) * first_dU

    # replace out of range dU values
    dpot_dU = np.where(dpot_flag == 'i', dpot_dU, first_dU_array)

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
parser.add_argument('dpot', type=argparse.FileType('w'))
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
    dpot_dU, dpot_flag = calculate_dpot(r, rdf_target_g, rdf_target_flag,
                                        rdf_current_g, rdf_current_flag,
                                        pot_current_U, pot_current_flag,
                                        args.kBT)

    # save dpot
    comment = "created by: {}".format(" ".join(sys.argv))
    saveto_table(args.dpot, r, dpot_dU, dpot_flag, comment)
