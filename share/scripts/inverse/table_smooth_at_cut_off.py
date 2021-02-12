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

import argparse
try:
    import numpy as np
except ImportError:
    print("Numpy is not installed, but needed for this script.")
from iie import find_nearest_ndx, saveto_table, readin_table

np.seterr(all='raise')


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
    ndx_co = find_nearest_ndx(r, cut_off)
    second_last_deriv = U[ndx_co-1] - U[ndx_co-2]
    shift = -1.0 * second_last_deriv - U[ndx_co-1]
    # modify up to second last value
    U_fixed[:ndx_co] += shift
    return U_fixed


def main():
    description = """\
    This script smoothes the potential close to the cut-off.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input', type=argparse.FileType('r'),
                        help='potential input')
    parser.add_argument('output', type=argparse.FileType('wb'),
                        help='potential output')
    parser.add_argument('--cut-off', type=float,
                        required=True, help='cut-off')
    args = parser.parse_args()
    # close writable files directly due to weird bug, where np.savetxt would
    # write empty file, use filename later.
    args.output.close()
    r, U, flag = readin_table(args.input)
    U_out = fix_U_near_cut_off_full(r, U, args.cut_off)
    saveto_table(args.output.name, r, U_out, flag, '')


if __name__ == '__main__':
    main()
