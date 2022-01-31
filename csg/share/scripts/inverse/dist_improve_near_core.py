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
try:
    import numpy as np
except ImportError:
    print("Numpy is not installed, but needed for this script")
    raise
from numpy.polynomial import Polynomial
if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.5+.")

np.seterr(all='raise')


def get_args(iie_args=None):
    """Define and parse command line arguments."""
    description = """
    Extrapolate an RDF close to the cores where its values are very small.

    It works by assuming a exponential or power form for the repulsive region of the PMF
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--verbose', dest='verbose',
                        help='save some intermeditary results',
                        action='store_const', const=True, default=False)
    parser.add_argument('--function', type=str, choices=['power', 'exponential'],
                        required=True,
                        help='Function used for PMF extrapolation.')
    parser.add_argument('--gmin', type=float, required=True,
                        help='lowest value to consider valid')
    parser.add_argument('--gmax', type=float, required=True,
                        help='highest value to consider valid')
    parser.add_argument('--in', type=argparse.FileType('r'),
                        required=True, dest='in_file',
                        help='RDF in file')
    parser.add_argument('--out', type=argparse.FileType('w'),
                        required=True, dest='out_file',
                        help='RDF out file')
    # parse
    args = parser.parse_args()
    return args


def process_input(args):
    """Process arguments and perform some checks."""
    # load input arrays
    r, g, flag = readin_table(args.in_file)
    # close writable file here because of np.savetxt bug
    args.out_file = args.out_file.name
    return r, g, flag


def improve_dist_near_core(r, g, pmf_function, fit_start_g, fit_end_g):
    g_extrap = g.copy()
    # find first peak
    # in some cases the second peak might be higher, then fit_start_g needs to
    # be smaller than valley
    g_max_ndx = np.argmax(g)
    # fit start: behind last point where g is smaller than fit_start_g
    fit_start_ndx = max(np.nonzero(g < fit_start_g)[0]) + 1
    # fit end: smallest index larger than fit_start_g, before first peak
    fit_end_ndx = min(np.nonzero(g[0:g_max_ndx] > fit_end_g)[0])
    if fit_end_ndx - fit_start_ndx < 3:
        raise Exception("less then three points found for fitting. This function needs "
                        f"a finer RDF in order to work. {fit_start_ndx} {fit_end_ndx}")
    # fitting
    if pmf_function == 'exponential':
        # fit r: ln(-ln(g)) with r: m * r + n
        # linearize
        def transform_x(x): return x
        def transform_y(y): return np.log(-np.log(y))
        def transform_y_inv(y): return np.exp(-np.exp(y))
    elif pmf_function == 'power':
        # fit ln(r): ln(-ln(g)) with m * r + n
        # linearize
        def transform_x(x): return np.log(x)
        def transform_y(y): return np.log(-np.log(y))
        def transform_y_inv(y): return np.exp(-np.exp(y))
    else:
        raise NotImplementedError("The value '{}' for pmf_function is not supported."
                                  "".format(pmf_function))
    # prepare data
    data_x = transform_x(r[fit_start_ndx:fit_end_ndx])
    data_y = transform_y(g[fit_start_ndx:fit_end_ndx])
    # fit
    fit = Polynomial.fit(data_x, data_y, 1)
    # use fit to extrap
    with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
        g_extrap[0:fit_start_ndx] = transform_y_inv(
            fit(transform_x(r[0:fit_start_ndx])))
    return g_extrap


def main():
    # get command line arguments
    args = get_args()
    # process and prepare input
    r, g, g_flag = process_input(args)
    # extrapolate
    g_improved = improve_dist_near_core(r, g, args.function, args.gmin,
                                        args.gmax)
    # modify flag?
    # save
    saveto_table(args.out_file, r, g_improved, g_flag, 'improved RDF')


if __name__ == '__main__':
    main()
