#!/usr/bin/env python3
"""Convert IMC matrix dS/du to a dg/du Jacobian."""
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
#
# Symbols:
# g: RDF
# U: potential
# dU: potential update (U_{k+1} - U_k)
#
# suffixes:
# _cur: current (of step k if currently doing iteration k)
# _tgt: target
# _mat: matrix (bead, bead) in the last two dimensions
# _vec: _mat matrix with all interactions flattened (vectorized)
# _2D: 2D matrix (which can be transformed to 4D)
# _flat: flat version of a vector, corresponds to matrix with _2D
#
# prefixes:
# ndx_: index

import argparse
import itertools
import numpy as np
import sys
import xml.etree.ElementTree as ET

from csg_functions import (
    calc_grid_spacing,
    calc_slices,
    gen_beadtype_property_array,
    gen_interaction_matrix,
    get_bead_types,
    get_density_dict,
    get_n_intra_dict,
    get_non_bonded,
    if_verbose_dump_io,
    make_matrix_2D,
    make_matrix_4D,
    read_all_tables,
    triangular_number,
    vectorize,
)

if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.5+.")

# raise all numpy errors. If errors are expected, use np.errstate()
np.seterr(all="raise")


def main():
    # get command line arguments
    args = get_args()
    # process and prepare input
    input_arrays, settings, imc_matrices, imc_indices = process_input(args)
    # convert and save
    convert_and_save_imc_matrix(input_arrays, settings, imc_matrices, imc_indices)


def get_args(local_args=None):
    """Define and parse command line arguments.

    If local_args is given, parse them instead of cmdlineargs.
    """
    description = "Calculate dg/du from dS/du."
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
        "--imc-matrix",
        type=argparse.FileType("r"),
        required=True,
        help="Text file with the IMC matrix",
    ),
    parser.add_argument(
        "--imc-index",
        type=argparse.FileType("r"),
        required=True,
        help="Text file with the IMC indices",
    ),
    parser.add_argument(
        "--volume",
        type=float,
        required=True,
        metavar="VOL",
        help="the volume of the box",
    )
    parser.add_argument(
        "--topol",
        type=argparse.FileType("r"),
        required=True,
        metavar="TOPOL",
        help="XML topology file",
    )
    parser.add_argument(
        "--options",
        type=argparse.FileType("r"),
        required=True,
        metavar="SETTINGS",
        help="XML settings file",
    )
    parser.add_argument(
        "--out",
        type=argparse.FileType("w"),
        required=True,
        help="Output .npz file with the Jacobian",
    )
    parser.add_argument(
        "--g-tgt-ext",
        type=str,
        required=True,
        metavar="RDF_TGT_EXT",
        help="extension of RDF target files",
    )
    parser.add_argument(
        "--g-cur-ext",
        type=str,
        required=True,
        metavar="RDF_CUR_EXT",
        help="extension of current RDF files",
    )
    parser.add_argument(
        "--improve-jacobian-onset",
        help="Change Jacobian slightly for better results in the RDF onset region",
        action="store_const",
        const=True,
        default=False,
    )
    parser.add_argument(
        "--onset-threshold",
        type=float,
        help="minimum value of g_tgt or g_cur up to which the IMC jacobian is improved",
    )
    parser.add_argument(
        "--cut-res",
        type=float,
        default=None,
        help="cut-off radius for each block",
    )
    # parse
    if local_args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(local_args)
    return args


def process_input(args):
    """Process arguments and perform some checks."""
    # args.options.read() can be called only once
    options = ET.fromstring(args.options.read())
    topology = ET.fromstring(args.topol.read())
    density_dict = get_density_dict(topology, args.volume)
    n_intra_dict = get_n_intra_dict(topology)

    # get non_bonded_dict
    non_bonded_dict = {nb_name: nb_ts for nb_name, nb_ts in get_non_bonded(options)}
    non_bonded_dict_inv = {v: k for k, v in non_bonded_dict.items()}
    if len(non_bonded_dict) != len(non_bonded_dict_inv):
        raise Exception(
            "Some non-bonded name was not unique or some non-bonded "
            "interactions had the same bead types."
        )

    # dict of table extensions
    table_infos = {
        "g_tgt": {
            "extension": args.g_tgt_ext,
            "check-grid": True,
            "assume-zero": [],
        },
        "g_cur": {
            "extension": args.g_cur_ext,
            "check-grid": True,
            "assume-zero": [],
        },
    }
    state_names_temp = ["."]  # read from current dir, dict will be collapsed later
    # read in all tables based on table_infos
    input_arrays, r0_removed, r = read_all_tables(
        state_names_temp,
        table_infos,
        non_bonded_dict,
    )
    # input_arrays structure has one level less if not multistate
    input_arrays = input_arrays["."]  # collapse
    del state_names_temp
    # quick access to r
    input_arrays["r"] = r
    # density and n_intra
    rhos = gen_beadtype_property_array(density_dict, non_bonded_dict)
    n_intra = gen_beadtype_property_array(n_intra_dict, non_bonded_dict)
    # settings
    # copy some settings directly from args
    args_to_copy = (
        "out",
        "volume",
        "verbose",
        "improve_jacobian_onset",
        "cut_res",
        "onset_threshold",
    )
    settings = {key: vars(args)[key] for key in args_to_copy}  # all mandatory
    settings["kBT"] = float(options.find("./inverse/kBT").text)
    settings["non-bonded-dict"] = non_bonded_dict
    settings["rhos"] = rhos
    settings["n_intra"] = n_intra
    settings["r0-removed"] = r0_removed

    if settings["improve_jacobian_onset"] and settings["onset_threshold"] is None:
        raise Exception(
            "If --improve-jacobian-onset is used, --onset-threshold has to be provided"
        )

    # load IMC matrix and index
    try:
        # close file(s), because we use np.load on file name
        args.imc_matrix.close()
        args.imc_index.close()
        imc_matrix = np.loadtxt(args.imc_matrix.name)
        imc_index = np.loadtxt(args.imc_index.name, dtype=str, ndmin=2)
    except (FileNotFoundError, ValueError):
        raise Exception("Can not load imc matrix/index file that was provided")

    return input_arrays, settings, imc_matrix, imc_index


@if_verbose_dump_io
def convert_and_save_imc_matrix(
    input_arrays, settings, imc_matrix, imc_index, verbose=False
):
    """Convert and save imc matrix to dg/du matrix and save as .npz.

    Args:
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        imc_matrices: IMC matrices to be converted
        verbose: save parameters and return of this and contained functions as numpy
                 file
    """
    # obtain r, Δr
    r = input_arrays["r"]
    Delta_r = calc_grid_spacing(r)
    # number of atom types
    n_t = len(settings["rhos"])
    # number of interactions excluding redundand ones (in(distinguishable))
    n_in = triangular_number(n_t)
    # number of grid points in r
    n_r = len(r)
    # connection non-bonded name <-> bead types
    non_bonded_dict = settings["non-bonded-dict"]
    # non_bonded_dict_inv = {v: k for k, v in settings["non-bonded-dict"].items()}
    # for each state calculate new matrix
    # index magic
    indices_imc = {}
    indices_new = {}
    # IMC matrix might contain r=0 values, remove them
    r0offset = 1 if settings["r0-removed"] else 0
    for i, (interaction, start_end) in enumerate(imc_index):
        start, end = map(int, start_end.split(":"))
        start_imc = start - 1  # votca counts from 1, Python counts from 0
        start_imc += r0offset  # remove leading r=0
        end_imc = end
        start_new = start - 1  # votca counts from 1, Python counts from 0
        start_new -= r0offset * i  # less values due to removed r=0
        end_new = end - r0offset * (i + 1)  # less values due to removed r=0
        indices_imc[interaction] = (start_imc, end_imc)
        indices_new[interaction] = (start_new, end_new)

    dudg_2D = np.zeros([n_r * n_in] * 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        for (nb1, index_imc1), (nb2, index_imc2) in itertools.product(
            indices_imc.items(), indices_imc.items()
        ):
            # Some lines to get the number of beads for interaction nb1
            # which relates to S in dS/du, because thate is wherere N1*N2
            # is the prefactor needed to obtain g from S
            beadtypes = tuple(non_bonded_dict[nb1])
            if len(beadtypes) == 1:
                bead_type1 = bead_type2 = beadtypes[0]
            else:
                bead_type1, bead_type2 = beadtypes
            all_beadtypes = get_bead_types(non_bonded_dict)
            bead_index1 = all_beadtypes.index(bead_type1)
            bead_index2 = all_beadtypes.index(bead_type2)
            N1 = (settings["rhos"] * settings["volume"])[bead_index1]
            N2 = (settings["rhos"] * settings["volume"])[bead_index2]

            # new indices
            index_new1 = indices_new[nb1]
            index_new2 = indices_new[nb2]

            # Kronecker delta for interactions with equal beads
            delta_alpha_beta = 1 if bead_index1 == bead_index2 else 0

            # dS/du -> dg/du
            dudg_2D[index_new1[0] : index_new1[1], index_new2[0] : index_new2[1]] = (
                1
                / settings["kBT"]
                * (delta_alpha_beta + 1)
                * settings["volume"]
                / (N1 * N2 * 4 * np.pi * Delta_r)
                * np.diag(1 / r**2)
                @ imc_matrix[
                    index_imc1[0] : index_imc1[1],
                    index_imc2[0] : index_imc2[1],
                ]
            )
    # make jacobian 4D
    dudg = make_matrix_4D(dudg_2D, n_r, n_r, n_in, n_in)

    # improve jacobian
    if settings["improve_jacobian_onset"]:
        dudg = improve_jacobian_onset(dudg, input_arrays, settings, verbose=verbose)

    # optionally cut Jacobian
    if settings["cut_res"]:
        cut, _ = calc_slices(r, settings["cut_res"])
        dudg = dudg[cut, cut, :, :]

    # save jacobian
    np.savez_compressed(settings["out"].name, jacobian=dudg)


@if_verbose_dump_io
def improve_jacobian_onset(jac_mat, input_arrays, settings, verbose=False):
    """Change Jacobian diagonal such that update becomes more stable.
    In the onset region the diagonal is approx. -k_B T g(r).
    This functions subtracts -k_B T g(r) (target or current) and adds
    -k_B T (g_cur + g_tgt)/2.

    Args:
        jac_mat: Jacobian to be improved (4D)
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file

    Returns:
        The improved Jacobian (4D)
    """
    r = input_arrays["r"]
    # n_r = len(r)
    # number of atom types
    n_t = len(settings["rhos"])
    # number of interactions including redundand ones
    n_i = triangular_number(n_t)

    # generate matrices and vectorize
    g_tgt_mat = gen_interaction_matrix(
        r, input_arrays["g_tgt"], settings["non-bonded-dict"]
    )
    g_tgt_vec = vectorize(g_tgt_mat)
    g_cur_mat = gen_interaction_matrix(
        r, input_arrays["g_cur"], settings["non-bonded-dict"]
    )
    g_cur_vec = vectorize(g_cur_mat)

    # copy jacobian
    jac_mat_improved = jac_mat.copy()

    # improve jacobian
    # loop over interactions
    for i in range(n_i):
        # onset region is up to where g_cur or g_tgt are larger then onset threshold
        onset_end = max(
            (g_cur_vec[:, i] > settings["onset_threshold"]).nonzero()[0][0],
            (g_tgt_vec[:, i] > settings["onset_threshold"]).nonzero()[0][0],
        )
        onset_region = slice(0, onset_end)

        # make all elements in the onset region zero
        jac_mat_improved[:, onset_region, :, i] = 0
        jac_mat_improved[onset_region, :, i, :] = 0

        # improve diagonal in onset retgion
        onset_diagonal = np.zeros(len(r[onset_region]))
        onset_diagonal = (
            -1
            / settings["kBT"]
            * ((g_cur_vec[onset_region, i] + g_tgt_vec[onset_region, i]) / 2)
        )

        # improve stability by converting zeros to small float
        onset_diagonal[onset_diagonal == 0.0] = -1e-37

        # change diagonal
        np.fill_diagonal(
            jac_mat_improved[onset_region, onset_region, i, i], onset_diagonal
        )
    # new diagonal
    return jac_mat_improved


if __name__ == "__main__":
    main()
