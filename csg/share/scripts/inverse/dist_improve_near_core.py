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
    readin_table,
    saveto_table,
)
import numpy as np
import sys

if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.5+.")

np.seterr(all="raise")


def main():
    # get command line arguments
    args = get_args()
    # process and prepare input
    r, g, g_flag = process_input(args)
    # extrapolate
    g_improved = improve_dist_near_core(r, g, args.gmin, args.gmax)
    # modify flag?
    # save
    saveto_table(args.out_file, r, g_improved, g_flag, "improved RDF")


def get_args(iie_args=None):
    """Define and parse command line arguments."""
    description = """
    Extrapolate an RDF close to the cores where its values are very small.

    It works by assuming a exponential + constant form for the repulsive region of the
    PMF.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        help="save some intermeditary results",
        action="store_const",
        const=True,
        default=False,
    )
    parser.add_argument(
        "--gmin",
        type=float,
        required=True,
        help="lowest value to consider valid",
    )
    parser.add_argument(
        "--gmax",
        type=float,
        required=True,
        help="highest value to consider valid",
    )
    parser.add_argument(
        "--in",
        type=argparse.FileType("r"),
        required=True,
        dest="in_file",
        help="RDF in file",
    )
    parser.add_argument(
        "--out",
        type=argparse.FileType("w"),
        required=True,
        dest="out_file",
        help="RDF out file",
    )
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


def gen_start_end_ndx(g, fit_start_g, fit_end_g):
    # find first peak
    # in some cases the second peak might be higher, then fit_start_g needs to
    # be smaller than valley
    g_max_ndx = np.argmax(g)
    # fit start: behind last point where g is smaller than fit_start_g
    fit_start_ndx = max(np.nonzero(g < fit_start_g)[0]) + 1
    # fit end: smallest index larger than fit_start_g, before first peak
    fit_end_ndx = min(np.nonzero(g[0:g_max_ndx] > fit_end_g)[0])
    return fit_start_ndx, fit_end_ndx


def f_exp_plus_const(x, a, b, c):
    return a * np.exp(-b * x) + c


def fit_exp_plus_const(x, y):
    """Fiting x, y with y = a * exp(-b * x) + c.

    To do this without scipy, we first fit the derivative and then find c
    in a second step.
    """
    Delta_x = x[1] - x[0]
    # derivative only depends on a and b
    # dy/dx = -ab exp(-b x)
    dydx = np.diff(y) / np.diff(x)
    # grid for derivative
    x_offset = x[:-1] + Delta_x / 2
    # log(-dy/dx) = log(ab) - bx
    try:
        log_m_dydx = np.log(-dydx)
    except FloatingPointError:
        raise Exception(
            """
The fitting of the PMF with a * exp(-b*x) + c is not possible, because it
is not monotonically decreasing.
Try better sampled RDFs or change the fitting range.
        """
        )
    # fit
    m_b, log_ab = np.polyfit(x=x_offset, y=log_m_dydx, deg=1)
    # unravel
    b = -m_b
    if b < 0:
        raise Exception(
            """
The fitting of the PMF with a * exp(-b*x) + c lead to a negative b.
This is non-physical. Likely your data in the fitting region is too noisy.
Try better sampled RDFs or change the fitting range.
        """
        )
    a = np.exp(log_ab) / b
    c = np.mean(y - f_exp_plus_const(x, a, b, 0))
    return a, b, c


def improve_dist_near_core(r, g, fit_start_g, fit_end_g):
    g_extrap = g.copy()
    fit_start_ndx, fit_end_ndx = gen_start_end_ndx(g, fit_start_g, fit_end_g)
    if fit_end_ndx - fit_start_ndx < 3:
        raise Exception(
            "less then three points found for fitting. This function needs "
            "a finer RDF or different fit limits in order to work. Current "
            f"fit indices: {fit_start_ndx} {fit_end_ndx}"
        )
    # prepare data
    data_x = r[fit_start_ndx:fit_end_ndx]
    data_y = -np.log(g[fit_start_ndx:fit_end_ndx])
    # fit
    abc = fit_exp_plus_const(data_x, data_y)
    # use fit to extrap
    with np.errstate(divide="ignore", invalid="ignore", under="ignore"):
        g_extrap[0:fit_start_ndx] = np.exp(-f_exp_plus_const(r[0:fit_start_ndx], *abc))
    return g_extrap


if __name__ == "__main__":
    main()
