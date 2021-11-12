#!/usr/bin/env python3
"""Extrapolate an RDF close to the core where the values are low and the sampling is
bad."""
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

import argparse
from csg_functions import (
    readin_table, saveto_table,
)
import sys
import numpy as np
if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.5+.")

np.seterr(all='raise')


def get_args(iie_args=None):
    """Define and parse command line arguments."""
    description = """
    Fix Δu, the potential update, near the cut-off where there are sometimes jumps with
    the Gauss-Newton method.

    It works by assuming constant force between the last three points of the potential.
    The last 'i' flag is assumed to be at the cut-off
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('in_file', type=argparse.FileType('r'), help='dpot in file')
    parser.add_argument('out_file', type=argparse.FileType('w'), help='dpot out file')
    # parse
    args = parser.parse_args()
    return args


def process_input(args):
    """Process arguments and perform some checks."""
    # load input arrays
    r, du, flag = readin_table(args.in_file)
    # close writable file here because of np.savetxt bug
    args.out_file = args.out_file.name
    return r, du, flag


def fix_dpot_near_cut_off(r, du, flag):
    du_shifted = du.copy()
    # index cut-off
    ndx_co = (flag == 'i').nonzero()[0][-1]
    # shift all but the last value by (value of last point + slope before that point)
    du_shifted[:ndx_co] -= du[ndx_co-1] + (du[ndx_co-1] - du[ndx_co-2])
    return du_shifted


def main():
    # get command line arguments
    args = get_args()
    # process and prepare input
    r, du, flag = process_input(args)
    # extrapolate
    du_fixed = fix_dpot_near_cut_off(r, du, flag)
    # save
    saveto_table(args.out_file, r, du_fixed, flag, 'fixed at cut-off')


if __name__ == '__main__':
    main()
