#!/usr/bin/env python3
"""Make a table flat behind the given cut-off."""
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
from csg_functions import readin_table, saveto_table, find_after_cut_off_ndx
import sys
import numpy as np

if not sys.version_info >= (3, 2):
    raise Exception("This script needs Python 3.2+.")

np.seterr(all="raise")


def get_args(iie_args=None):
    """Define and parse command line arguments."""
    description = """
    Make the y column of the potential zero behind the cut-off.

    The rest of the values are shifted, such that the value at the cut-off is zero.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "cut_off", type=float, help="cut-off after which will be leveled"
    )
    parser.add_argument("in_file", type=argparse.FileType("r"), help="dpot in file")
    parser.add_argument("out_file", type=argparse.FileType("w"), help="dpot out file")
    # parse
    args = parser.parse_args()
    return args


def process_input(args):
    """Process arguments and perform some checks."""
    # load input arrays
    r, u, flag = readin_table(args.in_file)
    # close writable file here because of np.savetxt bug
    args.out_file.close()
    # continue with filename
    args.out_file = args.out_file.name
    cut_off = args.cut_off
    return r, u, flag, cut_off


def level_table_after_cut_off(r, u, cut_off):
    # index cut-off
    ndx_aco = find_after_cut_off_ndx(r, cut_off)
    # do notghing if cut_off outside of table
    if ndx_aco >= len(r):
        return u
    # copy array
    u_leveled = u.copy()
    # level after cut_off
    u_leveled[ndx_aco:] = 0.0
    # shift up to cut_off
    u_leveled[:ndx_aco] -= u[ndx_aco - 1]
    return u_leveled


def main():
    # get command line arguments
    args = get_args()
    # process and prepare input
    r, u, flag, cut_off = process_input(args)
    # extrapolate
    u_leveled = level_table_after_cut_off(r, u, cut_off)
    # save
    saveto_table(args.out_file, r, u_leveled, flag, "fixed at cut-off")


if __name__ == "__main__":
    main()
