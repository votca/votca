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
    gen_beadtype_property_array,
    gen_interaction_matrix,
    get_density_dict,
    get_n_intra_dict,
    get_non_bonded,
    if_verbose_dump_io,
    make_matrix_2D,
    make_matrix_4D,
    read_all_tables,
    vectorize,
)

if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.5+.")

# constants
BAR_PER_MD_PRESSURE = 16.6053904  # md pressure is kJ/mol/nm³ as in Gromacs
F_COULOMB = 138.935458  # electric conversion factor: V = f*q_1*q_2/r

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
        nargs="+",
        default=None,
        help="Text file with the IMC matrix. Multiple files for multistate.",
    ),
    parser.add_argument(
        "--imc-index",
        type=argparse.FileType("r"),
        nargs="+",
        default=None,
        help="Text file with the IMC indices. Multiple files for multistate.",
    ),
    parser.add_argument(
        "--volume",
        type=float,
        required=True,
        nargs="+",
        metavar="VOL",
        help="the volume of the box. Multiple if multistate",
    )
    parser.add_argument(
        "--topol",
        type=argparse.FileType("r"),
        required=True,
        nargs="+",
        metavar="TOPOL",
        help="XML topology file. Multiple if multistate",
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
        nargs="+",
        help="Output .npz file with the Jacobian. Multiple files for multistate.",
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
        "--g-cur-raw-ext",
        type=str,
        default=None,
        metavar="RDF_CUR_RAW_EXT",
        help="extension of current raw RDF files",
    )
    parser.add_argument(
        "--g-tgt-raw-ext",
        type=str,
        default=None,
        metavar="RDF_CUR_RAW_EXT",
        help="extension of current raw RDF files",
    )
    parser.add_argument(
        "--improve-jacobian-onset",
        help="Change Jacobian slightly for better results in the RDF onset region",
        action="store_const",
        const=True,
        default=False,
    )
    parser.add_argument(
        "--multistate",
        dest="multistate",
        action="store_const",
        const=True,
        default=False,
        help="enable multistate method",
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
    # multistate settings
    multistate = args.multistate
    if multistate:
        state_names = options.find("./inverse/multistate/state_names").text.split()

    # get topology, density_dict, and n_intra_dict
    if multistate:
        topology = [ET.fromstring(top_file.read()) for top_file in args.topol]
        density_dict = [
            get_density_dict(top, vol) for top, vol in zip(topology, args.volume)
        ]
        n_intra_dict = [get_n_intra_dict(top) for top in topology]  # prob. indep.
    else:
        topology = ET.fromstring(args.topol[0].read())
        density_dict = get_density_dict(topology, args.volume[0])
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
    if args.g_cur_raw_ext is not None:
        table_infos = {
            **table_infos,
            "g_cur_raw": {
                "extension": args.g_cur_raw_ext,
                "check-grid": True,
                "assume-zero": [],
            },
        }
    if args.g_tgt_raw_ext is not None:
        table_infos = {
            **table_infos,
            "g_tgt_raw": {
                "extension": args.g_tgt_raw_ext,
                "check-grid": True,
                "assume-zero": [],
            },
        }
    if multistate:
        state_names_temp = state_names
    else:
        state_names_temp = ["."]  # read from current dir, dict will be collapsed later
    # read in all tables based on table_infos
    input_arrays, r0_removed, r = read_all_tables(
        state_names_temp,
        table_infos,
        non_bonded_dict,
    )
    # input_arrays structure has one level less if not multistate
    if not multistate:
        assert len(input_arrays) == 1
        input_arrays = input_arrays["."]  # collapse
    del state_names_temp
    # quick access to r
    input_arrays["r"] = r
    # density and n_intra
    if multistate:
        rhos = [gen_beadtype_property_array(dd, non_bonded_dict) for dd in density_dict]
        n_intra = [
            gen_beadtype_property_array(nd, non_bonded_dict) for nd in n_intra_dict
        ]
    else:
        rhos = gen_beadtype_property_array(density_dict, non_bonded_dict)
        n_intra = gen_beadtype_property_array(n_intra_dict, non_bonded_dict)
    # settings
    # copy some settings directly from args
    args_to_copy = (
        "out",
        "volume",
        "verbose",
        "improve_jacobian_onset",
    )
    settings = {key: vars(args)[key] for key in args_to_copy}  # all mandatory
    settings["kBT"] = float(options.find("./inverse/kBT").text)
    settings["non-bonded-dict"] = non_bonded_dict
    settings["rhos"] = rhos
    settings["n_intra"] = n_intra
    settings["r0-removed"] = r0_removed

    # load Joacbian and index
    imc_matrices = []  # one per state
    imc_indices = []
    for imc_matrix, imc_index in zip(args.imc_matrix, args.imc_index):
        try:
            # close file(s), because we use np.load on file name
            imc_matrix.close()
            imc_index.close()
            imc_matrices.append(np.loadtxt(imc_matrix.name))
            imc_indices.append(np.loadtxt(imc_index.name, dtype=str, ndmin=2))
        except (FileNotFoundError, ValueError):
            raise Exception("Can not load imc matrix file(s) that were provided")
    if multistate:
        assert len(imc_matrices) == len(state_names)
    else:
        assert len(imc_matrices) == 1

    # multistate settings
    settings["multistate"] = multistate
    if multistate:
        settings["state_names"] = state_names
        settings["state_weights"] = list(
            map(float, options.find("./inverse/multistate/state_weights").text.split())
        )
        settings["state_kBTs"] = list(
            map(float, options.find("./inverse/multistate/state_kBTs").text.split())
        )

    return input_arrays, settings, imc_matrices, imc_indices


@if_verbose_dump_io
def convert_and_save_imc_matrix(
    input_arrays, settings, imc_matrices, imc_indices, verbose=False
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
    # number of interactions including redundand ones
    n_i = int(n_t**2)
    # number of grid points in r
    n_r = len(r)
    # for each state calculate new matrix
    for imc_matrix, imc_index, outfile in zip(
        imc_matrices, imc_indices, settings["out"]
    ):
        indices = {}
        for i, (interaction, start_end) in enumerate(imc_index):
            start, end = map(int, start_end.split(":"))
            start -= 1  # votca counts from 1, I count from 0
            end += 1  # using end for slicing
            indices[interaction] = (start, end, i)

        # IMC matrix might contain 0, remove it
        r0offset = 1 if settings["r0-removed"] else 0
        dudg_2D = np.zeros([n_r * n_i] * 2)
        with np.errstate(divide="ignore", invalid="ignore"):
            for (nb1, index1), (nb2, index2) in itertools.product(
                indices.items(), indices.items()
            ):
                N1 = (settings["rhos"] * settings["volume"])[index1[2]]
                N2 = (settings["rhos"] * settings["volume"])[index2[2]]
                # dS/du -> dg/du
                # indices make sure that r=0 is excluded
                dudg_2D[index1[0] : index1[1], index2[0] : index2[1]] = (
                    1
                    / settings["kBT"]
                    / (N1 * N2 / settings["volume"] / 2 * 4 * np.pi * Delta_r)
                    * np.diag(1 / r**2)
                    @ imc_matrix[
                        index1[0] + r0offset * (index1[2] + 1) : index1[1],
                        index2[0] + r0offset * (index1[2] + 1) : index2[1],
                    ]
                )
        # make jacobian 4D
        dudg = make_matrix_4D(dudg_2D, n_r, n_r, n_i, n_i)

        # improve jacobian
        if settings["improve_jacobian_onset"]:
            dudg = improve_jacobian_onset(dudg, input_arrays, settings, verbose=verbose)

        # save jacobian
        np.savez_compressed(outfile.name, jacobian=dudg)


@if_verbose_dump_io
def improve_jacobian_onset(J, input_arrays, settings, verbose=False):
    """Change Jacobian diagonal such that update becomes more stable.
    In the onset region the diagonal is approx. -k_B T g(r).
    This functions subtracts -k_B T g(r) (target or current) and adds
    -k_B T (g_cur + g_tgt)/2.

    Args:
        dudg: Jacobian to be improved
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file

    Returns:
        The improved Jacobian
    """
    r = input_arrays["r"]
    n_r = len(r)
    # number of atom types
    n_t = len(settings["rhos"])
    # number of interactions including redundand ones
    n_i = int(n_t**2)
    # generate matrices and vectorize
    g_tgt_mat = gen_interaction_matrix(
        r, input_arrays["g_tgt"], settings["non-bonded-dict"]
    )
    g_tgt_vec = vectorize(g_tgt_mat)
    g_cur_mat = gen_interaction_matrix(
        r, input_arrays["g_cur"], settings["non-bonded-dict"]
    )
    g_cur_vec = vectorize(g_cur_mat)
    # when dist_improve was used, subtract raw
    if "g_cur_raw" in input_arrays.keys():
        vec_to_subract = vectorize(
            gen_interaction_matrix(
                r, input_arrays["g_cur_raw"], settings["non-bonded-dict"]
            )
        ).T.flatten()
    elif "g_tgt_raw" in input_arrays.keys():
        vec_to_subract = vectorize(
            gen_interaction_matrix(
                r, input_arrays["g_tgt_raw"], settings["non-bonded-dict"]
            )
        ).T.flatten()
    else:
        vec_to_subract = g_cur_vec.T.flatten()
    J_improved_2D = make_matrix_2D(J.copy())
    # TODO: IMC does not distinguish u_ab, u_ba, _vec does. Implement and use hvec?
    J_diag = np.diag(J_improved_2D)
    # subtract -g_cur then add average of (-g_cur, -g_tgt) for better stability
    # TODO: This would need another version, when use_target_matrix is used
    J_diag_improved = J_diag + 1 / settings["kBT"] * (
        +vec_to_subract - (g_cur_vec.T.flatten() + g_tgt_vec.T.flatten()) / 2
    )
    np.fill_diagonal(J_improved_2D, J_diag_improved)
    J_improved = make_matrix_4D(J_improved_2D, n_r, n_r, n_i, n_i)
    return J_improved


if __name__ == "__main__":
    main()
