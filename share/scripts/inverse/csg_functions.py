#!/usr/bin/env python3
"""Multi purpose script for Iterative Integral Equation methods."""
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

from collections import defaultdict
import sys
try:
    import numpy as np
except ImportError:
    print("Numpy is not installed, but needed for the the CSG python functions")
    raise
if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.5+.")


def readin_table(filename):
    """Read in votca table."""
    table_dtype = {'names': ('x', 'y', 'y_flag'),
                   'formats': ('f', 'f', 'U2')}
    x, y, y_flag = np.loadtxt(filename, dtype=table_dtype, comments=['#', '@'],
                              unpack=True)
    return x, y, y_flag


def saveto_table(filename, x, y, y_flag, comment=""):
    """Save votca table."""
    data = np.zeros((len(x),), dtype='f, f, U2')
    data['f0'] = x
    data['f1'] = y
    data['f2'] = y_flag
    np.savetxt(filename, data, header=comment, fmt=['%e', '%e', '%s'])


def calc_grid_spacing(grid, relative_tolerance=0.01):
    """Returns the spacing of an equidistant 1D grid."""
    diffs = np.diff(grid)
    if abs((max(diffs) - min(diffs)) / max(diffs)) > relative_tolerance:
        raise Exception('the grid is not equidistant')
    return np.mean(diffs)


def test_calc_grid_spacing():
    """Check spacing of some grid."""
    grid = np.linspace(0, 2 * np.pi, num=361)
    grid_spacing = calc_grid_spacing(grid)
    assert np.allclose(grid_spacing, np.pi/180)


def fourier(r, f):
    """Compute the radially 3D FT of a radially symmetric function.

    The frequency grid is also returned.  Some special extrapolations are used
    to make the results consistent. This function is isometric meaning it can
    be used to calculate the FT and the inverse FT.  That means inputs can also
    be k and f_hat which results in r and f.

    Args:
        r: Input grid. Must be evenly spaced. Can start at zero or at Δr, but nowhere
            else.
        f: Input function. Must have same length as r and correspond to its values.

    Returns:
        (k, f_hat): The reciprocal grid and the FT of f.

    """
    Delta_r = calc_grid_spacing(r)
    r0_added = False
    if np.isclose(r[0], Delta_r):
        r = np.concatenate(([0], r))
        f = np.concatenate(([0], f))
        r0_added = True
    elif np.isclose(r[0], 0.0):
        pass
    else:
        raise Exception('this function can not handle this input')
    # if the input is even, np.fft.rfftfreq would end with the Nyquist frequency.
    # But there the imaginary part of the FT is always zero, so we alwas append a zero
    # to obtain a odd grid.
    if len(r) % 2 == 0:  # even
        r = np.concatenate((r, [r[-1]+Delta_r]))
        f = np.concatenate((f, [0]))
        n = (len(r)-1)*2-1
    else:  # odd
        n = len(r)*2-1
    k = np.fft.rfftfreq(n=n, d=Delta_r)
    with np.errstate(divide='ignore', invalid='ignore'):
        f_hat = -2 / k / 1 * Delta_r * np.imag(np.fft.rfft(r * f, n=n))
    if r0_added:
        f_hat = f_hat[1:]
        k = k[1:]
    return k, f_hat


def test_fourier():
    """Check that Fourier function is invertible."""
    r = np.linspace(1, 100, 100)
    f = np.random.random(100)
    k, f_hat = fourier(r, f)
    r_, f_ = fourier(k, f_hat)
    assert np.allclose(r, r_)
    assert np.allclose(f, f_)


def fourier_all(r, table_mat):
    """Fourier along first dimension of table matrix (all interactions)."""
    k, _ = fourier(r, table_mat[:, 0, 0])  # get length
    table_hat_mat = np.empty((len(k), *table_mat.shape[1:]))
    table_hat_mat.fill(np.nan)
    for b1 in range(table_mat.shape[1]):
        for b2 in range(table_mat.shape[2]):
            y = table_mat[:, b1, b2]
            _, y_hat = fourier(r, y)
            table_hat_mat[:, b1, b2] = y_hat
    return k, table_hat_mat


def gen_fourier_matrix(r, fourier_function):
    """Make a fourier matrix."""
    fourier_matrix = np.identity(len(r))
    for col_index, col in enumerate(fourier_matrix.T):
        _, fourier_matrix.T[col_index] = fourier_function(r, col)
    return fourier_matrix


def find_nearest_ndx(array, value):
    """Find index of array where closest to value."""
    array = np.asarray(array)
    ndx = (np.abs(array - value)).argmin()
    return ndx


def test_find_nearest_ndx():
    """Check finding the nearest index."""
    tests = [
        ([0, 1, 2, 3], -1, 0),
        ([0, 1, 2, 3], 4.9, 3),
        ([0, 1, 2, 3], 1.49, 1),
        ([0, 1, 2, 3], 1.51, 2),
    ]
    for grid, val, ndx in tests:
        assert find_nearest_ndx(grid, val) == ndx


def find_after_cut_off_ndx(array, cut_off):
    """
    Find index of array after given cut_off.

    Assumes array is sorted. Used for finding first index after cut_off.

    """
    array = np.asarray(array)
    ndx_closest = find_nearest_ndx(array, cut_off)
    if np.isclose(array[ndx_closest], cut_off):
        return ndx_closest + 1
    if array[-1] < cut_off:
        return len(array)
    ndx = np.where(array > cut_off)[0][0]
    return ndx


def test_find_after_cut_off_ndx():
    """Check finding the index after a value."""
    tests = [
        ([0, 1, 2, 3], 1.0, 2),
        ([0, 1, 2, 3], 1.01, 2),
        ([0, 1, 2, 3], 1.99, 2),
        ([0, 1, 2, 3], 3.5, 4),
    ]
    for grid, val, ndx in tests:
        print(find_after_cut_off_ndx(grid, val), ndx)
        assert find_after_cut_off_ndx(grid, val) == ndx


def r0_removal(*arrays):
    """
    Remove the first element from a list of arrays.

    Only does so if the first array starts with 0.

    """
    r0_removed = False
    if np.isclose(arrays[0][0], 0.0):
        r0_removed = True
        arrays = tuple(map(lambda a: a[1:], arrays))
    return r0_removed, arrays


def get_non_bonded(options_xml):
    """Yield tuple (name, bead_types) for each non-bonded interaction.

    bead_types is a set of the interaction's bead types."""
    for non_bonded in options_xml.findall('non-bonded'):
        type_set = frozenset({non_bonded.find('type1').text,
                              non_bonded.find('type2').text})
        yield non_bonded.find('name').text, type_set


def get_density_dict(topol_xml, volume):
    """Return the densities of all bead types as a dict."""
    density_dict = defaultdict(lambda: 0.0)
    for molecule in topol_xml.find('molecules').findall('molecule'):
        for bead in molecule.findall('bead'):
            density_dict[bead.attrib['type']] += int(molecule.attrib['nmols']) / volume
    density_dict = dict(density_dict)
    return density_dict


def get_n_intra_dict(topol_xml):
    """Return the number beads per molecules."""
    n_intra_dict = defaultdict(lambda: 0)
    for molecule in topol_xml.find('molecules').findall('molecule'):
        for bead in molecule.findall('bead'):
            n_intra_dict[bead.attrib['type']] += 1
    n_intra = dict(n_intra_dict)
    return n_intra


def get_bead_types(non_bonded_dict):
    """Return a sorted list of bead types."""
    bead_types = {bead
                  for value in non_bonded_dict.values()
                  for bead in value}
    bead_types = sorted(list(bead_types))
    return bead_types


def gen_interaction_matrix(r, interaction_dict, non_bonded_dict):
    bead_types = get_bead_types(non_bonded_dict)
    non_bonded_dict_inv = {v: k for k, v in non_bonded_dict.items()}
    interaction_matrix = np.empty((len(r), len(bead_types), len(bead_types)))
    interaction_matrix.fill(np.nan)
    for b1, bead1 in enumerate(bead_types):
        for b2, bead2 in enumerate(bead_types):
            interaction_name = non_bonded_dict_inv[frozenset({bead1, bead2})]
            interaction_matrix[:, b1, b2] = interaction_dict[interaction_name]['y']
    return interaction_matrix


# inverse of the above
def gen_interaction_dict(r, interaction_matrix, non_bonded_dict):
    bead_types = get_bead_types(non_bonded_dict)
    non_bonded_dict_inv = {v: k for k, v in non_bonded_dict.items()}
    interaction_dict = {}
    for b1, bead1 in enumerate(bead_types):
        for b2, bead2 in enumerate(bead_types):
            interaction_name = non_bonded_dict_inv[frozenset({bead1, bead2})]
            interaction_dict[interaction_name] = {'x': r,
                                                  'y': interaction_matrix[:, b1, b2]}
    return interaction_dict


def gen_beadtype_property_array(property_dict, non_bonded_dict):
    bead_types = get_bead_types(non_bonded_dict)
    try:
        property_array = np.array([property_dict[bt] for bt in bead_types])
    except KeyError:
        raise Exception("Could not construct density array. Inconsistency between "
                        "topology and options file?")
    return property_array


def gauss_newton_constrained(A, C, b, d):
    """Do a gauss-newton update, but eliminate Cx=d first."""
    m, n = A.shape
    p, n_ = C.shape
    assert n == n_
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
        print("WARNING: solution of Gauss-Newton update determined fully "
              "by constraints.")
        x_elim = []
    else:
        x_elim = np.linalg.solve(A_elim.T @ A_elim, A_elim.T @ b_elim)
    if p == 0:
        # no constraints
        x = x_elim
    else:
        x_pivot = (d[i] - np.delete(C, pivot, 1) @ x_elim) / C[i, pivot]
        x = np.insert(x_elim, pivot, x_pivot)
    return x


def test_gauss_newton_constrained():
    """Check Gauss-Newton with some simple cases."""
    tests = [
        (np.identity(10), np.ones((1, 10)), np.ones(10), np.array(2), [0.2]*10),
        (np.identity(5), np.array([[0, 0, 1, 0, 0]]), np.ones(5), np.array(2),
         [1, 1, 2, 1, 1]),
        (np.array([[1, 0], [1, 1]]), np.zeros((0, 2)), np.ones(2), np.array([]),
         [1.0, 0.0]),
        (np.array([[1, 0], [1, 1]]), np.array([[0, 1]]), np.ones(2), np.array(0.1),
         [0.95, 0.1]),
    ]
    for A, C, b, d, x in tests:
        assert np.allclose(x, gauss_newton_constrained(A, C, b, d))


def upd_flag_g_smaller_g_min(flag, g, g_min):
    """
    Update the flag to 'o' for small RDF.

    Take a flag list, copy it, and set the flag to 'o'utside if g is smaller
    g_min.

    """
    flag_new = flag.copy()
    for i, gg in enumerate(g):
        if gg < g_min:
            flag_new[i] = 'o'
    return flag_new


def upd_flag_by_other_flag(flag, other_flag):
    """
    Update a flag array by another flag array.

    Take a flag list, copy it, and set the flag to 'o'utside where some
    other flag list is 'o'.

    """
    flag_new = flag.copy()
    for i, of in enumerate(other_flag):
        if of == 'o':
            flag_new[i] = 'o'
    return flag_new


def gen_flag_isfinite(U):
    """
    Generate a flag list based on if the elements of U are finite.
    """
    return np.where(np.isfinite(U), ['i'] * len(U), ['o'] * len(U))


def extrapolate_dU_left_constant(dU, dU_flag):
    """
    Extrapolate the potential update in the core region by a constant value.

    The first valid value, determined by the flag, is used
    """
    dU_extrap = dU.copy()
    # find first valid dU value
    first_dU_index = np.where(dU_flag == 'i')[0][0]
    first_dU = dU[first_dU_index]

    # replace out of range dU values with constant first value
    dU_extrap = np.where(dU_flag == 'i', dU, first_dU)
    return dU_extrap
