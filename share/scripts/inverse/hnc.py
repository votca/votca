#!/usr/bin/env python3
""" fooo """
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
# Symbols:
# g: RDF
# U: potential
# dU: potential update (U_{k+1} - U_k)
#
# suffixes:
# _cur: current or of step k if currently doing iteration k
# _tgt: target
# _ce: core_end (where RDF gets > 0)
# _co: cut_off
# _sc: single_component
#
# prefixes:
# ndx_: index


import argparse
import sys
import numpy as np


BAR_PER_MD_PRESSURE = 16.6053904
G_MIN = 1e-10
np.seterr(all='raise')


def readin_table(filename):
    """read in votca table"""
    table_dtype = {'names': ('x', 'y', 'y_flag'),
                   'formats': ('f', 'f', 'S1')}
    x, y, y_flag = np.loadtxt(filename, dtype=table_dtype, comments=['#', '@'],
                              unpack=True)
    return x, y, y_flag


def saveto_table(filename, x, y, y_flag, comment=""):
    """save votca table"""
    data = np.column_stack((x.T, y.T, y_flag.T))
    np.savetxt(filename, data, header=comment, fmt='%s')


def compare_grids(grid_a, grid_b):
    """check two grids for no point differing more than 0.1 picometer"""
    if np.any(grid_a - grid_b > 0.0001):
        raise Exception("Different grids!")


def calc_grid_spacing(grid, relative_tolerance=0.01):
    """returns the spacing of an equidistant 1D grid.
    fails if not equidistant grid."""
    diffs = np.diff(grid)
    if abs((max(diffs) - min(diffs)) / max(diffs)) > relative_tolerance:
        raise Exception('the grid is not equidistant')
    return np.mean(diffs)


def fourier(r, y, omega):
    """fourier transform real data y on grid r into reciprocal space omega"""
    Delta_r = calc_grid_spacing(r)
    y_hat = np.zeros_like(omega)
    with np.errstate(divide='ignore', invalid='ignore', under='ignore'):
        for i, omega_i in enumerate(omega):
            y_hat[i] = (2 / omega_i * Delta_r
                        * np.sum(r * y * np.sin(2 * np.pi * omega_i * r)))
    return y_hat


def gen_fourier_matrix(r_nz, omega_nz, fourier_function):
    """make a fourier matrix"""
    fourier_matrix = np.identity(len(omega_nz))
    for col_index, col in enumerate(fourier_matrix.T):
        fourier_matrix.T[col_index] = fourier_function(r_nz, col, omega_nz)
    return fourier_matrix


def find_nearest_ndx(array, value):
    """find index of array where closest to value"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def gauss_newton_constrained(A, C, b, d):
    """do a gauss-newton update, but eliminate Cx=d first"""
    m, n = A.shape
    p, n = C.shape
    b.shape = (m)
    d.shape = (p)

    if p > 1:
        raise Exception("not implemented for p > 1")

    A_elim = A.copy()
    b_elim = b.copy()
    for i in range(p):
        pivot = np.argmax(abs(C[i]))  # find max value of C
        A_elim = A - (np.ones_like(A) * A[:, pivot][:, np.newaxis]
                      * C[i] / C[i, pivot])
        b_elim = b - A[:, pivot] * d[i] / C[i, pivot]
        A_elim = np.delete(A_elim, pivot, 1)
    if p == n:
        print("Warning: solution determined fully by constraints")
        x_elim = []
    else:
        x_elim = np.linalg.solve(np.matmul(A_elim.T, A_elim),
                                 np.matmul(A_elim.T, b_elim))
    if p == 0:
        # no constraints
        x = x_elim
    else:
        x_pivot = (d[i] - np.matmul(np.delete(C, pivot, 1),
                                    x_elim)) / C[i, pivot]
        x = np.insert(x_elim, pivot, x_pivot)
    return x


def calc_dU_hncgn_sc(r, g_cur, g_tgt, kBT, density,
                     g_min, cut_off, constraints,
                     verbose):
    """Apply the (constrained) HNCGN method."""

    # for convenience grids (real and reciprocal) are without the zero value
    r, g, g_tgt = r[1:], g_cur[1:], g_tgt[1:]
    # there are different regions in r used in the method
    #              |       crucial     |                     # regions
    # |   core     |                 nocore               |
    #  ---------core_end------------cut_off-----------r[-1]  # distances
    # 0----------ndx_ce--------------ndx_co----------len(r)  # indices
    #
    # Δ from the paper equals nocore
    # Δ' from the paper equals crucial
    # note: Vector w is in region crucial, but with one element less, because
    # of the antiderivative operator A0
    ndx_ce = np.where(g_tgt > g_min)[0][0]
    ndx_co = find_nearest_ndx(r, cut_off)
    crucial = slice(ndx_ce, ndx_co+1)
    nocore = slice(ndx_ce, len(r))
    Delta_r = calc_grid_spacing(r)
    if verbose:
        print("core_end:", ndx_ce, r[ndx_ce])
        print("cut_off:", ndx_co, r[ndx_co])
        print("r_end:", len(r)-1, r[-1])
        print("r", 0, len(r), min(r), max(r))
        print("crucial:", ndx_ce, ndx_co+1, min(r[crucial]), max(r[crucial]))
        print("nocore:", ndx_ce, len(r), min(r[nocore]), max(r[nocore]))

    # reciprocal grid up to Nyquist frequency
    nyquist_freq = 1 / Delta_r / 2
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
    T = np.matmul(np.linalg.inv(F), np.matmul(H, F))
    # D matrix
    D = np.diag(g_tgt[ndx_ce:])
    # U matrix
    U = -kBT * np.linalg.inv(D) + kBT * T[ndx_ce:, ndx_ce:]
    # A0 matrix
    A0 = Delta_r * np.triu(np.ones((len(r[nocore]), len(r[crucial])-1)), k=0)
    # Jacobian matrix
    J = np.matmul(np.linalg.inv(U), A0)
    # constraint matrix and vector
    C = np.zeros((len(constraints), len(r[crucial])-1))
    d = np.zeros(len(constraints))
    # build constraint matrix and vector from constraints
    for c, constraint in enumerate(constraints):
        if constraint['type'] == 'pressure':
            # current pressure
            p = constraint['current'] / BAR_PER_MD_PRESSURE
            # target pressure
            p_tgt = constraint['target'] / BAR_PER_MD_PRESSURE
            # g_tgt(r_{i+1})
            g_tgt_ip1 = g_tgt[crucial][1:]
            # g_tgt(r_{i})
            g_tgt_i = g_tgt[crucial][:-1]
            # r_{i+1}
            r_ip1 = r[crucial][1:]
            # r_{i}
            r_i = r[crucial][:-1]
            # l vector
            ll = (g_tgt_i + g_tgt_ip1) * (r_ip1**4 - r_i**4)
            ll *= 1/12 * np.pi * rho**2
            # set C row and d element
            C[c, :] = ll
            d[c] = p_tgt - p
        else:
            raise Exception("not implemented constraint type")
    # residuum vector
    res = g_tgt - g
    # switching to notation of Gander et al. for solving
    A = J
    b = res[nocore]
    w = gauss_newton_constrained(A, C, b, d)
    # dU
    dU = np.matmul(A0, w)
    # fill core with nans
    dU = np.concatenate((np.full(ndx_ce + 1, np.nan), dU))
    # dump files
    if verbose:
        np.savez_compressed('hncgn-arrays.npz', A=A, b=b, C=C, d=d,
                            H=H, T=T, D=D, U=U, A0=A0, J=J)
    return dU


def extrapolate_U_constant(dU, dU_flag):
    """Extrapolate the potential in the core region by a constant value.

    Returns dU. U_{k+1} = U_k + dU is done by Votca."""
    dU_extrap = dU.copy()
    # find first valid dU value
    first_dU_index = np.where(dU_flag == 'i')[0][0]
    first_dU = dU[first_dU_index]

    # replace out of range dU values with constant first value
    dU_extrap = np.where(dU_flag == 'i', dU, first_dU)
    return dU_extrap


def extrapolate_U_power(r, dU, U, g_tgt, g_min, kBT, verbose):
    """Extrapolate the potential in the core region including
    the point, where the RDF becomes larger than g_min or where
    the new potential is convex. A power function is used.
    The PMF is fitted, not U+dU. The fit is then shifted such
    that the graph is monotonous.

    Extrapolation is done, because p-HNCGN often has artifacs at
    the core, especially when pressure is far off.

    Returns dU. U_{k+1} = U_k + dU is done by Votca."""
    # make copy
    dU_extrap = dU.copy()
    # calc PMF
    with np.errstate(divide='ignore', over='ignore'):
        pmf = np.nan_to_num(-kBT * np.log(g_tgt))
    # index first minimum
    ndx_fm = np.where(np.nan_to_num(np.diff(pmf)) > 0)[0][0]
    # index core end
    ndx_ce = np.where(g_tgt > g_min)[0][0]
    # fit pmf region
    fit_region = slice(ndx_ce, ndx_ce + 3)
    # fit pmf with power function a*x^b
    pmf_shift = -pmf[ndx_fm] + 0.01
    fit_x = np.log(r[fit_region])
    fit_y = np.log(pmf[fit_region] + pmf_shift)
    b, log_a = np.polyfit(fit_x, fit_y, 1)
    a = np.exp(log_a)
    if verbose:
        print('pmf fit a*x^b. Coefficients a, b:', a, b)
    # r without zero at start
    with np.errstate(divide='ignore', over='ignore'):
        pmf_fit = np.nan_to_num(a * r**b - pmf_shift)

    # region to extrapolate
    ndx_ex1 = ndx_ce + 1
    ndx_ex2 = np.where(np.nan_to_num(np.diff(np.diff(U + dU))) > 0)[0][0]
    ndx_ex = max(ndx_ex1, ndx_ex2)
    # extrapolate
    U_extrap = U + dU
    U_extrap[:ndx_ex] = pmf_fit[:ndx_ex] + (U_extrap[ndx_ex] - pmf_fit[ndx_ex])
    dU_extrap = U_extrap - U
    return dU_extrap


def shift_U_cutoff_zero(dU, r, U, cut_off):
    """Make potential zero at and beyond cut-off"""
    dU_shift = dU.copy()
    # shift dU to be zero at cut_off and beyond
    ndx_co = find_nearest_ndx(r, cut_off)
    U_cut_off = U[ndx_co] + dU[ndx_co]
    dU_shift -= U_cut_off
    dU_shift[ndx_co:] = -1 * U[ndx_co:]
    return dU_shift


def fix_U_near_cut_off_full(dU, r, U, cut_off):
    """Modify the potential close to the cut-off in
    a way, such that it is more smooth. The derivative
    of the potential between the last two points will
    be equal to the derivative between the two points
    before. The original last two points of dU are
    therefore ignored.

    This also helps agains an artifact of p-HNCGN,
    where the last value of dU is a spike."""
    dU_fixed = dU.copy()
    ndx_co = find_nearest_ndx(r, cut_off)
    U_new = U + dU
    second_last_deriv = U_new[ndx_co-1] - U_new[ndx_co-2]
    shift = -1.0 * second_last_deriv - U_new[ndx_co-1]
    # modify up to second last value
    dU_fixed[:ndx_co] += shift
    return dU_fixed


def calc_dU_ihnc_sym_mol(r, g_tgt, g_cur, G_minus_g, n, kBT, rho,
                         hncn_variant=False):
    """calculates an update step dU for the potential from g_tgt, g_cur, and
    G_minus_g  using IHNC (or HNCN) step for molecules with n identical beads.
    This means it works for two-bead hexane, but is expected to fail for
    three-bead hexane"""
    # reciprocal space with radial symmetry ω
    omega = np.arange(1, len(r)) / (2 * max(r))
    # difference of rdf to target 'f'
    g_prime = g_tgt - g_cur
    g_prime_hat = fourier(r, g_prime, omega)
    # pair correlation function 'h'
    h_hat = fourier(r, g_cur - 1, omega)
    # intramolecular distribution G - g
    G_minus_g_hat = fourier(r, G_minus_g, omega)
    # auxiliary variables a, b, e
    a_hat = rho * G_minus_g_hat
    # reduces to 1 for a_hat = 0 or n = 1
    b_hat = 1 + a_hat * (2 - 2/n) + a_hat**2 * (1/n * (n - 2) + 1/n**2)
    # reduces to rho for a_hat = 0 or n = 1
    e_hat = rho * (a_hat - a_hat/n + 1)
    # dc/dg g' = c'
    # derived in sympy
    c_prime_hat = b_hat * g_prime_hat / (b_hat + e_hat * h_hat)**2
    c_prime = fourier(omega, c_prime_hat, r)
    # dU from HNC
    with np.errstate(divide='ignore', invalid='ignore'):
        if hncn_variant:
            dU = kBT * (-(g_prime / g_tgt) + g_prime - c_prime)
        else:
            dU = kBT * (np.log(g_cur / g_tgt) + g_prime - c_prime)
    return dU


def calc_U_hnc_sym_mol(r, g_cur, G_minus_g, n, kBT, rho):
    """calculates U from g and G_minus_g for molecules with n identical beads.
    This means it works for two-bead hexane, but is expected to fail for
    three-bead hexane"""
    # reciprocal space with radial symmetry ω
    omega = np.arange(1, len(r)) / (2 * max(r))
    # total correlation function h
    h = g_cur - 1
    h_hat = fourier(r, h, omega)
    # intramolecular distribution G - g
    G_minus_g_hat = fourier(r, G_minus_g, omega)
    # auxiliary variables a, b, e
    a_hat = rho * G_minus_g_hat
    # reduces to 1 for a_hat = 0 or n = 1
    b_hat = 1 + a_hat * (2 - 2/n) + a_hat**2 * (1/n * (n - 2) + 1/n**2)
    # reduces to rho for a_hat = 0 or n = 1
    e_hat = rho * (a_hat - a_hat/n + 1)
    # direct correlation function c from OZ including intramolecular
    # interactions
    c_hat = h_hat / (b_hat + e_hat * h_hat)
    c = fourier(omega, c_hat, r)
    # U from HNC
    with np.errstate(divide='ignore', invalid='ignore'):
        U = kBT * (-np.log(g_cur) + h - c)
    return U


def upd_flag_g_smaller_g_min(flag, g, g_min):
    """Take a flag list, copy it, and set the flag to 'o'utside if g is smaller
    g_min"""
    flag_new = flag.copy()
    for i, gg in enumerate(g):
        if gg < g_min:
            flag_new[i] = 'o'
    return flag_new


def upd_flag_by_other_flag(flag, other_flag):
    """Take a flag list, copy it, and set the flag to 'o'utside where some
    other flag list is 'o'"""
    flag_new = flag.copy()
    for i, of in enumerate(other_flag):
        if of == 'o':
            flag_new[i] = 'o'
    return flag_new


def upd_U_const_first_flag_i(U, flag):
    """Take a potential list, copy it, and set the potential at places where
    flag is 'o' to the first value where flag is 'i'"""
    U_new = U.copy()
    # find first valid U value
    first_U_index = np.where(flag == 'i')[0][0]
    first_U = U_new[first_U_index]
    # replace out of range U values
    U_new = np.where(flag == 'i', U_new, first_U)
    return U_new


def upd_U_zero_beyond_cut_off(r, U, cut_off):
    """Take a potential list, copy it, and shift U to be zero at cut_off and
    beyond"""
    U_new = U.copy()
    index_cut_off = np.searchsorted(r, cut_off)
    U_cut_off = U_new[index_cut_off]
    U_new -= U_cut_off
    U_new[index_cut_off:] = 0
    return U_new


def main():
    # always raise error

    description = """\
    This script calculatess U or ΔU with the HNC methods.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-v', '--verbose', dest='verbose',
                        help='save some intermeditary results',
                        action='store_const', const=True, default=False)
    # subparsers
    subparsers = parser.add_subparsers(dest='subcommand', required=True)
    parser_invert_hnc = subparsers.add_parser('invert_hnc',
                                              help='potential guess using HNC')
    parser_ihnc = subparsers.add_parser('ihnc',
                                        help='potential update IHNC')
    parser_hncgn = subparsers.add_parser('hncgn',
                                         help='potential update HNCGN')

    # all subparsers
    for pars in [parser_invert_hnc, parser_ihnc, parser_hncgn]:
        pars.add_argument('--g-tgt', type=argparse.FileType('r'),
                          nargs='+', required=True,
                          metavar=('X-X.dist.tgt', 'X-Y.dist.tgt'),
                          help='RDF target files')
        pars.add_argument('--kBT', type=float, required=True, help='')
        pars.add_argument('--densities', type=float,
                          nargs='+', required=True,
                          metavar=('rho_X', 'rho_Y'),
                          help='list of number densities of beads')
        pars.add_argument('--n-intra', type=int,
                          nargs='+', required=True,
                          metavar=('n_X', 'n_Y'),
                          help='number of beads per molecule')
        # todo: this might need to be split up for multicomponent
        pars.add_argument('--cut-off', type=float, required=True,
                          help='cutt-off (co) of potential')
        pars.add_argument('--U-out', type=argparse.FileType('wb'),
                          nargs='+', required=True,
                          metavar=('X-X.dpot.new', 'X-Y.dpot.new'),
                          help='U or ΔU output files')
    # intial potential guess subparsers
    for pars in [parser_invert_hnc]:
        pars.add_argument('--G-tgt', type=argparse.FileType('r'),
                          nargs='+', required=False,
                          metavar=('X-X-incl.dist.tgt', 'X-Y-incl.dist.tgt'),
                          help='RDF (including intramolecular) target files')
    # update potential subparsers
    for pars in [parser_ihnc, parser_hncgn]:
        pars.add_argument('--g-cur', type=argparse.FileType('r'),
                          nargs='+', required=True,
                          metavar=('X-X.dist.cur', 'X-Y.dist.cur'),
                          help='RDF current files')
        pars.add_argument('--G-cur', type=argparse.FileType('r'),
                          nargs='*', required=False,
                          metavar=('X-X-incl.dist.cur', 'X-Y-incl.dist.cur'),
                          help='RDF (including intramolecular) current files')
        pars.add_argument('--U-cur', type=argparse.FileType('r'),
                          nargs='+', required=True,
                          metavar=('X-X.pot.cur', 'X-Y.pot.cur'),
                          help='potential current files')
    # HNCGN only options
    parser_hncgn.add_argument('--pressure-constraint',
                              dest='pressure_constraint',
                              type=str, default=None)
    parser_hncgn.add_argument('--extrap-near-core',
                              dest='extrap_near_core',
                              type=str, choices=['none', 'constant', 'power'])
    parser_hncgn.add_argument('--fix-near-cut-off',
                              dest='fix_near_cut_off',
                              type=str, choices=['none', 'full-deriv'])

    args = parser.parse_args()

    # close writable files directly due to weird bug, where np.savetxt would
    # write empty file, use filename later.
    # I should report it, but on the other side, argparse opened files do not
    # get closed at any point, so this is better
    for f in args.U_out:
        f.close()
    args.U_out = [f.name for f in args.U_out]

    # infering variables from input
    n_beads = len(args.densities)
    # nr. of elements in triangular matrix incl. diagonal
    n_interactions = (n_beads * (n_beads + 1)) // 2

    # some checks on input
    # test for same number of interactions
    file_arguments_names = ['g_tgt', 'G_tgt', 'g_cur', 'G_cur', 'U_cur',
                            'U_out']
    file_arguments = [(argname, vars(args)[argname])
                      for argname in file_arguments_names
                      if vars(args).get(argname) is not None]
    file_arguments_wrong = [(argname, flist)
                            for argname, flist in file_arguments
                            if len(flist) != n_interactions]
    for argname, flist in file_arguments_wrong:
        raise Exception(f"""You provided N = {n_beads} densities, therefore
                        there should be (N * (N + 1)) // 2 = {n_interactions}
                        files for {argname}, but {[f.name for f in flist]} was
                        provided""")
    # multicomponent not implemented
    if any(len(files) != 1 for files in [args.g_tgt, args.densities]):
        raise Exception('not implemented for multiple components!')
    if len(args.n_intra) > 1:
        raise Exception('not implemented for multiple components!')
    # if n_intra > 1, also needs G
    # NOT TRUE: may just want to initiate without G
    """
    G_arguments_names = ['G_tgt', 'G_cur']
    if args.n_intra[0] > 1 and not any(vars(args).get(G_name)
                                       for G_name in G_arguments_names):
        raise Exception('If n_intra is larger than 1, you should also '
                        'provide some G')
    """
    # todo: if n_intra == 1, check if G close to g

    # todo for multicomponent: check order of input and output by filename
    # todo for multicomponent: allow not existing X-Y? particles would overlap
    # todo for multicomponent: do not allow same bead on different
    # moleculetypes

    # load input arrays
    input_arrays = {}  # input_arrays['g_tgt'][0]['x']
    input_arguments_names = ['g_tgt', 'G_tgt', 'g_cur', 'G_cur', 'U_cur']
    input_arguments = [(argname, vars(args)[argname])
                       for argname in input_arguments_names
                       if vars(args).get(argname) is not None]
    for argname, flist in input_arguments:
        input_arrays[argname] = []
        for i, f in enumerate(flist):
            x, y, flag = readin_table(f)
            input_arrays[argname].append({'x': x, 'y': y, 'flag': flag})

    # todo: compare grids of all
    r = input_arrays['g_tgt'][0]['x']
    # compare_grids(r, input_arrays['g_cur'][0]['x'])
    # compare grids of G_cur and G_tgt with g_tgt in smart way

    # calculate G_minus_g
    if 'G_tgt' in input_arrays:
        g_name = 'g_tgt'
        G_name = 'G_tgt'
    elif 'G_cur' in input_arrays:
        g_name = 'g_cur'
        G_name = 'G_cur'
    else:
        print('no intramolecular correlations provided, assuming there are '
              'none')
        # this will result in zero arrays in input_arrays['G_minus_g'] in the
        # below loop
        g_name = 'g_tgt'
        G_name = 'g_tgt'
    input_arrays['G_minus_g'] = []
    for g_dict, G_dict in zip(input_arrays[g_name], input_arrays[G_name]):
        G_minus_g = np.zeros_like(g_dict['y'])
        G_minus_g[:G_dict['y'].shape[0]] = (G_dict['y']
                                            - g_dict['y']
                                            [:G_dict['y'].shape[0]])
        input_arrays['G_minus_g'].append({'y': G_minus_g})

    # invert_hnc
    if args.subcommand in ['invert_hnc']:
        U1 = calc_U_hnc_sym_mol(r,
                                input_arrays['g_tgt'][0]['y'],
                                input_arrays['G_minus_g'][0]['y'],
                                args.n_intra[0], args.kBT, args.densities[0])
        U_flag1 = np.array(['i'] * len(r))
        U_flag2 = upd_flag_g_smaller_g_min(U_flag1,
                                           input_arrays['g_tgt'][0]['y'],
                                           G_MIN)
        U2 = upd_U_const_first_flag_i(U1, U_flag2)
        U3 = upd_U_zero_beyond_cut_off(r, U2, args.cut_off)
        comment = "created by: {}".format(" ".join(sys.argv))
        saveto_table(args.U_out[0], r, U3, U_flag2, comment)

    # ihnc
    if args.subcommand in ['ihnc']:
        dU1 = calc_dU_ihnc_sym_mol(r,
                                   input_arrays['g_tgt'][0]['y'],
                                   input_arrays['g_cur'][0]['y'],
                                   input_arrays['G_minus_g'][0]['y'],
                                   args.n_intra[0], args.kBT,
                                   args.densities[0])
        dU_flag1 = np.array(['i'] * len(r))
        dU_flag2 = upd_flag_g_smaller_g_min(dU_flag1,
                                            input_arrays['g_tgt'][0]['y'],
                                            G_MIN)
        dU_flag3 = upd_flag_g_smaller_g_min(dU_flag2,
                                            input_arrays['g_cur'][0]['y'],
                                            G_MIN)
        dU_flag4 = upd_flag_by_other_flag(dU_flag3,
                                          input_arrays['U_cur'][0]['flag'])

        dU2 = upd_U_const_first_flag_i(dU1, dU_flag4)
        U_temp = upd_U_zero_beyond_cut_off(r,
                                           dU2 + input_arrays['U_cur'][0]['y'],
                                           args.cut_off)
        dU3 = U_temp - input_arrays['U_cur'][0]['y']
        comment = "created by: {}".format(" ".join(sys.argv))
        saveto_table(args.U_out[0], r, dU3, dU_flag4, comment)

    # hncgn
    if args.subcommand in ['hncgn']:
        # parse constraints
        constraints = []
        if args.pressure_constraint is not None:
            p_target = float(args.pressure_constraint.split(',')[0])
            p_current = float(args.pressure_constraint.split(',')[1])
            constraints.append({'type': 'pressure', 'target': p_target,
                                'current': p_current})

        # calc dU_pure
        dU_pure = calc_dU_hncgn_sc(r,
                                   input_arrays['g_cur'][0]['y'],
                                   input_arrays['g_tgt'][0]['y'],
                                   args.kBT,
                                   args.densities[0], G_MIN, args.cut_off,
                                   constraints, args.verbose)

        # set dU_flag to 'o' inside the core
        dU_flag = np.where(np.isnan(dU_pure), 'o', 'i')

        # select extrapolation
        if args.extrap_near_core == 'none':
            dU_extrap = np.nan_to_num(dU_pure)
        elif args.extrap_near_core == 'constant':
            dU_extrap = extrapolate_U_constant(dU_pure, dU_flag)
        elif args.extrap_near_core == 'power':
            dU_extrap = extrapolate_U_power(r, dU_pure,
                                            input_arrays['U_cur'][0]['y'],
                                            input_arrays['g_tgt'][0]['y'],
                                            G_MIN, args.kBT, args.verbose)
        else:
            raise Exception("unknown extrapolation scheme for inside and near "
                            "core region: " + args.extrap_near_core)
        # shifts to correct potential after cut-off
        dU_shift = shift_U_cutoff_zero(dU_extrap, r,
                                       input_arrays['U_cur'][0]['y'],
                                       args.cut_off)
        # shifts to correct potential near cut-off
        if args.fix_near_cut_off == 'none':
            dU = dU_shift.copy()
        elif args.fix_near_cut_off == 'full-deriv':
            dU = fix_U_near_cut_off_full(dU_shift, r,
                                         input_arrays['U_cur'][0]['y'],
                                         args.cut_off)
        else:
            raise Exception("unknown fix scheme for near cut-off: "
                            + args.fix_near_cut_off)

        if args.verbose:
            np.savez_compressed('hncgn-dU.npz', r=r, dU_pure=dU_pure,
                                dU_extrap=dU_extrap, dU_shift=dU_shift)
        comment = "created by: {}".format(" ".join(sys.argv))
        saveto_table(args.U_out[0], r, dU, dU_flag, comment)


if __name__ == '__main__':
    main()
