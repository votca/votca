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
#
# Symbols:
# g: RDF
# u: potential
#
# suffixes:
# _cur: current (of step k if currently doing iteration k)
# _tgt: target
# _mat: matrix (bead, bead) in the last two dimensions
# _vec: _mat matrix with all interactions flattened (vectorized), corresponds to matrix
#        with _mat
# _2D: 2D matrix (which can be transformed to 4D)
# _flat: flat version of a vector, corresponds to matrix with _2D
#
# prefixes:
# ndx_: index

import argparse
import itertools
import sys
import xml.etree.ElementTree as ET
import numpy as np

from imc_matrix_to_jacobian import improve_jacobian_onset
from csg_functions import (
    calc_slices,
    cut_matrix_inverse,
    extrapolate_Delta_u_left_constant,
    fourier,
    fourier_all,
    gen_beadtype_property_array,
    gen_flag_isfinite,
    gen_fourier_matrix,
    gen_interaction_dict,
    gen_interaction_matrix,
    get_charge_dict,
    get_density_dict,
    get_intra_needed,
    get_n_intra_dict,
    get_non_bonded,
    if_verbose_dump_io,
    kron_2D,
    make_matrix_2D,
    make_matrix_4D,
    read_all_tables,
    saveto_table,
    transpose,
    vectorize_full,
)

if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.5+.")

# constants
F_COULOMB = 138.935458  # electric conversion factor: V = f*q_1*q_2/r
# raise all numpy errors. If errors are expected, use np.errstate()
np.seterr(all="raise")


def main():
    # get command line arguments
    args = get_args()
    # process and prepare input
    input_arrays, settings = process_input(args)
    # guess potential from distribution
    if settings["subcommand"] == "potential_guess":
        output_arrays = potential_guess(
            input_arrays, settings, verbose=settings["verbose"]
        )
        # save all potentials to table files
        save_tables(output_arrays, settings)
    # calculate dc/dh
    elif settings["subcommand"] == "dcdh":
        if settings["out"] is not None:
            calc_and_save_dcdh(input_arrays, settings, verbose=settings["verbose"])
    # newton update
    elif settings["subcommand"] == "jacobian":
        calc_and_save_jacobian(input_arrays, settings, verbose=settings["verbose"])


def get_args(iie_args=None):
    """Define and parse command line arguments.

    If iie_args is given, parse them instead of cmdlineargs.
    """
    description = "Calculate u(r), dc/du, or dg/du  with Integral Equations."
    parser = argparse.ArgumentParser(description=description)
    # subparsers
    subparsers = parser.add_subparsers(dest="subcommand")
    parser_pot_guess = subparsers.add_parser(
        "potential_guess", help="potential guess from inverting integral equation"
    )
    parser_dcdh = subparsers.add_parser(
        "dcdh",
        help="calculate the dc/dh matrix",
    )
    parser_jacobian = subparsers.add_parser(
        "jacobian",
        help="calculate dh/du, the Jacobian",
    )
    # all subparsers
    for pars in (parser_pot_guess, parser_dcdh, parser_jacobian):
        pars.add_argument(
            "-v",
            "--verbose",
            dest="verbose",
            help="save some intermeditary results",
            action="store_const",
            const=True,
            default=False,
        )
        pars.add_argument(
            "--volume",
            type=float,
            required=True,
            metavar="VOL",
            help="The volume of the box.",
        )
        pars.add_argument(
            "--topol",
            type=argparse.FileType("r"),
            required=True,
            metavar="TOPOL",
            help="XML topology file.",
        )
        pars.add_argument(
            "--options",
            type=argparse.FileType("r"),
            required=True,
            metavar="SETTINGS",
            help="XML settings file",
        )
        pars.add_argument(
            "--g-tgt-ext",
            type=str,
            required=True,
            metavar="RDF_TGT_EXT",
            help="extension of RDF target files",
        )
        pars.add_argument(
            "--out",
            type=str,
            required=True,
            metavar="U_OUT_EXT",
            help="extension of u or Δu files or full filename of dcdh "
            "matrix. If 'none' there will be no output.",
        )
        pars.add_argument(
            "--g-tgt-intra-ext",
            type=str,
            metavar="RDF_TGT_INTRA_EXT",
            help="extension of intramol. RDF target files",
        )
        pars.add_argument(
            "--kBT",
            type=float,
            required=True,
            metavar="kBT",
            dest="kBT",
            help="Temperature times k_B",
        )
    # pot_guess and jacobian need closure
    for pars in (parser_pot_guess, parser_jacobian):
        pars.add_argument(
            "--closure",
            type=str,
            choices=["hnc", "py"],
            required=True,
            help="Closure equation to use for the OZ equation",
        )
    # matrix generation needs cut-residual
    for pars in (parser_dcdh, parser_jacobian):
        pars.add_argument(
            "--cut-residual",
            type=float,
            required=True,
            help="The cut-off for the residues; will determine size for jac/dcdh",
        )
    # Jacobian needs current RDF for diagnoal and optionally also intra RDF and tgt_dcdh
    parser_jacobian.add_argument(
        "--g-cur-ext",
        type=str,
        required=True,
        metavar="RDF_CUR_EXT",
        help="extension of current RDF files",
    )
    parser_jacobian.add_argument(
        "--g-cur-intra-ext",
        type=str,
        metavar="RDF_CUR_INTRA_EXT",
        help="extension of current intramol. RDF files",
    )
    parser_jacobian.add_argument(
        "--tgt-dcdh",
        type=argparse.FileType("r"),
        default=None,
        help=(
            ".npz file with dc/dh from target distributions. "
            "If provided, will be used. "
            "Otherwise the jacobian will be calculated from "
            "current distributions."
        ),
    )
    parser_jacobian.add_argument(
        "--improve-jacobian-onset",
        help="Change Jacobian slightly for better results in the RDF onset region",
        action="store_const",
        const=True,
        default=False,
    )
    parser_jacobian.add_argument(
        "--onset-thresholds",
        nargs=2,
        type=float,
        default=None,
        help="two RDF values that determine the range in the onset region where the "
        "jacobian is changed from the original to the approximation.",
    )
    # potential guess can remove Coulomb term
    parser_pot_guess.add_argument(
        "--subtract-coulomb",
        dest="subtract_coulomb",
        help="remove Coulomb term from potential guess",
        action="store_const",
        const=True,
        default=False,
    )
    # parse
    if iie_args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(iie_args)
    # check for subcommand
    if args.subcommand is None:
        parser.print_help()
        raise Exception("subcommand needed")
    return args


def process_input(args):
    """Process arguments and perform some checks."""
    # args.options.read() can be called only once
    options = ET.fromstring(args.options.read())
    args.options.close()
    # get topology, density_dict, and n_intra_dict
    topology = ET.fromstring(args.topol.read())
    args.topol.close()
    density_dict = get_density_dict(topology, args.volume)
    n_intra_dict = get_n_intra_dict(topology)
    # get charge_dict
    if args.subcommand == "potential_guess":
        charge_dict = get_charge_dict(topology)
    else:
        charge_dict = None
    # get needed intramolecular interactions
    # all topologies should contain all molecules
    intra_needed = get_intra_needed(topology)
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
        "g_tgt": {"extension": args.g_tgt_ext, "check-grid": True, "assume-zero": []}
    }
    # if potential guess or update and tgt_jacobian we need the target intramolecular
    # RDFs
    if args.subcommand in ("potential_guess", "dcdh"):
        table_infos = {
            **table_infos,
            "G_minus_g_tgt": {
                "extension": args.g_tgt_intra_ext,
                "check-grid": True,
                "assume-zero": set(
                    [
                        non_bonded_dict_inv[beadset]
                        for beadset in (set(non_bonded_dict.values()) - intra_needed)
                    ]
                ),
            },
        }
    # if jacobian, we need the current RDFs for the diagonal
    if args.subcommand == "jacobian":
        table_infos = {
            **table_infos,
            "g_cur": {
                "extension": args.g_cur_ext,
                "check-grid": True,
                "assume-zero": [],
            },
        }
        # if not target jacobian we need the current intramolecular RDFs
        if args.tgt_dcdh is None:
            table_infos = {
                **table_infos,
                "G_minus_g_cur": {
                    "extension": args.g_cur_intra_ext,
                    "check-grid": True,
                    "assume-zero": set(
                        [
                            non_bonded_dict_inv[beadset]
                            for beadset in (
                                set(non_bonded_dict.values()) - intra_needed
                            )
                        ]
                    ),
                },
            }
    # load input arrays
    input_arrays = {}  # will hold all input data
    # Structure: [table_name, non_bonded_name, xyflag]
    # Multistate structure: [state, table_name, non_bonded_name, xyflag]
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
    # process input further
    rhos = gen_beadtype_property_array(density_dict, non_bonded_dict)
    n_intra = gen_beadtype_property_array(n_intra_dict, non_bonded_dict)
    # settings
    # copy some settings directly from args
    args_to_copy = (
        "closure",
        "verbose",
        "subcommand",
        "subtract_coulomb",
        "out",
        "kBT",
        "improve_jacobian_onset",
        "onset_thresholds",
    )
    settings = {key: vars(args)[key] for key in args_to_copy if key in vars(args)}
    settings["non-bonded-dict"] = non_bonded_dict
    settings["rhos"] = rhos
    settings["n_intra"] = n_intra
    settings["r0-removed"] = r0_removed

    if args.subcommand == "jacobian":
        if settings["improve_jacobian_onset"] and settings["onset_thresholds"] is None:
            raise Exception(
                "If --improve-jacobian-onset is used, "
                "--onset-thresholds has to be provided."
            )

    # determine dc/dh buffer
    if args.subcommand == "jacobian":
        if args.tgt_dcdh is None:
            settings["tgt_dcdh"] = None
        else:
            # reuse dc/dh
            # close file(s), because we use np.load on file name
            args.tgt_dcdh.close()
            try:
                settings["tgt_dcdh"] = np.load(args.tgt_dcdh.name)["dcdh"]
            except (FileNotFoundError, ValueError):
                raise Exception("Can not load tgt_dcdh file that was provided")
    # determine cut-off xml path
    if args.subcommand == "potential_guess":
        cut_off_path = "./inverse/initial_guess/cut_off"
        try:
            settings["cut_off_pot"] = float(options.find(cut_off_path).text)
        except (AttributeError, ValueError):
            raise Exception(
                cut_off_path + " must be a float in settings.xml for integral "
                "equation methods"
            )
    elif args.subcommand in ("dcdh", "jacobian"):
        settings["cut_off_res"] = args.cut_residual
    # stuff for subtracting coulomb potential
    if args.subcommand == "potential_guess":
        settings["charge_dict"] = charge_dict
        settings["non_bonded_dict"] = non_bonded_dict

    return input_arrays, settings


@if_verbose_dump_io
def potential_guess(input_arrays, settings, verbose=False):
    """Calculate potential guess based on symmetry adapted RISM-OZ and closure.

    Args:
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file

    Returns:
        dictionary of potentials including flags to be saved
    """
    # obtain r
    r = input_arrays["r"]
    n_r = len(r)
    # cut-off
    cut, tail = calc_slices(r, settings["cut_off_pot"], verbose=False)
    n_c = len(r[cut])
    if n_c * 2 > n_r:
        print(
            "WARNING: max is smaller than twice of cut_off. This will lead to "
            "artifacts in the Jacobian."
        )
    # prepare matrices
    g_mat = gen_interaction_matrix(
        r, input_arrays["g_tgt"], settings["non-bonded-dict"]
    )
    h_mat = g_mat - 1
    k, h_hat_mat = fourier_all(r, h_mat)
    G_minus_g_mat = gen_interaction_matrix(
        r, input_arrays["G_minus_g_tgt"], settings["non-bonded-dict"]
    )
    _, G_minus_g_hat_mat = fourier_all(r, G_minus_g_mat)
    # perform actual math
    u_mat = calc_u_matrix(
        r,
        k,
        g_mat,
        h_hat_mat,
        G_minus_g_hat_mat,
        settings["rhos"],
        settings["n_intra"],
        settings["kBT"],
        settings["closure"],
        verbose=settings["verbose"],
    )
    # subtract Coulomb, extrapolate, and save potentials
    output_arrays = {}
    for non_bonded_name, u_dict in gen_interaction_dict(
        r, u_mat, settings["non-bonded-dict"]
    ).items():
        u = u_dict["y"]
        u_flag = gen_flag_isfinite(u)
        # subtract Coulomb
        if settings["subtract_coulomb"]:
            beads = tuple(settings["non_bonded_dict"][non_bonded_name])
            bead1, bead2 = (beads[0], beads[0]) if len(beads) == 1 else beads
            q1 = settings["charge_dict"][bead1]
            q2 = settings["charge_dict"][bead2]
            u_Coulomb = F_COULOMB * q1 * q2 / r
            u -= u_Coulomb
        # make tail zero. It is spoiled on the last half from inverting OZ.
        # careful: slices refer to arrays before reinserting r=0 values!
        u[cut] -= u[cut][-1]
        u[tail] = 0
        u_flag[tail] = "o"
        # add value at r=0 if it was removed
        if settings["r0-removed"]:
            r_out = np.concatenate(([0.0], r))
            u = np.concatenate(([np.nan], u))
            u_flag = np.concatenate((["o"], u_flag))
        else:
            r_out = r
        # change NaN in the core region to first valid value
        u = extrapolate_Delta_u_left_constant(u, u_flag)
        output_arrays[non_bonded_name] = {"x": r_out, "y": u, "flag": u_flag}
    return output_arrays


def save_tables(output_arrays, settings):
    """Save each entry in output_arrays to a table file."""
    comment = "created by: {}".format(" ".join(sys.argv))
    if settings["out"] == "none":
        return None
    for non_bonded_name, output_dict in output_arrays.items():
        fn = non_bonded_name + "." + settings["out"]
        saveto_table(
            fn, output_dict["x"], output_dict["y"], output_dict["flag"], comment
        )


@if_verbose_dump_io
def calc_u_matrix(
    r,
    k,
    g_mat,
    h_hat_mat,
    G_minus_g_hat_mat,
    rhos,
    n_intra,
    kBT,
    closure,
    verbose=False,
):
    """
    Calculate a potential u using integral equation theory.

    Args:
        r: Distance grid.
        g_mat: matrix of RDF
        h_hat_mat: matrix of Fourier of TCF
        G_minus_g_mat: matrix of Fourier of intramolecular RDF
        rhos: array of densities of the bead types
        n_intra: array with number of bead per molecule
        kBT: Boltzmann constant times temperature.
        closure: OZ-equation closure ('hnc' or 'py').
        verbose: output calc_u_matrix.npz

    Returns:
        matrix of the calculated potentias.
    """
    # calculate direct correlation function
    c_mat = calc_c_matrix(r, k, h_hat_mat, G_minus_g_hat_mat, rhos, n_intra, verbose)
    with np.errstate(divide="ignore", invalid="ignore"):
        if closure == "hnc":
            u_mat = kBT * (-np.log(g_mat) + (g_mat - 1) - c_mat)
        elif closure == "py":
            u_mat = kBT * np.log(1 - c_mat / g_mat)
        else:
            raise Exception("closure has to be hnc or py")
    return u_mat


@if_verbose_dump_io
def calc_c_matrix(r, k, h_hat_mat, G_minus_g_hat_mat, rhos, n_intra, verbose=False):
    """Calculate the direct correlation function c from g for all interactions."""
    # row sum of ω, after Bertagnolli and my own notes
    Omega_hat_mat = gen_Omega_hat_mat(G_minus_g_hat_mat, rhos, n_intra)
    # H_hat_mat after Bertagnolli
    H_hat_mat = adapt_reduced_matrix(h_hat_mat, n_intra)
    # Rho_mat after Bertagnolli
    Rhos = rhos / n_intra
    # intermediate step
    # have to transpose to solve x·a = b with numpy by solving a'·x' = b'
    H_over_Omega_plus_rho_H = transpose(
        np.linalg.solve(
            transpose(Omega_hat_mat + np.diag(Rhos) @ H_hat_mat), transpose(H_hat_mat)
        )
    )
    # direct correlation function C from symmetry reduced OZ
    C_hat_mat = np.linalg.solve(Omega_hat_mat, H_over_Omega_plus_rho_H)
    # c_hat from C_hat
    c_hat_mat = unadapt_reduced_matrix(C_hat_mat, n_intra)
    # c from c_hat
    _, c_mat = fourier_all(k, c_hat_mat)
    return c_mat


@if_verbose_dump_io
def calc_jacobian(input_arrays, settings, verbose=False):
    """
    Calculate the Jacobian using integral equations

    Args:
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file

    Returns:
        The Jacobian™ and its inverse
    """
    # obtain r
    r = input_arrays["r"]
    # number of atom types
    n_t = len(settings["rhos"])
    # number of interactions
    n_i = int(n_t**2)
    # number of grid points in r
    n_r = len(r)
    # slices
    cut, _ = calc_slices(r, settings["cut_off_res"], verbose=False)
    n_c = len(r[cut])
    # generate matrices and cut them directly
    g_tgt_mat = gen_interaction_matrix(
        r, input_arrays["g_tgt"], settings["non-bonded-dict"]
    )
    g_cur_mat = gen_interaction_matrix(
        r, input_arrays["g_cur"], settings["non-bonded-dict"]
    )
    # which distributions to use for dc/dh
    # using cur is the original Newton-Raphson root finding method
    # using tgt is a method similar to Newton's but with slope calculated at the root
    # the latter is taken from the input, is calculated at step_000 once
    if settings["tgt_dcdh"] is not None:
        dcdh = settings["tgt_dcdh"]
    else:
        # generate dc/dh, invert, cut it, and invert again
        G_minus_g_cur_mat = gen_interaction_matrix(
            r, input_arrays["G_minus_g_cur"], settings["non-bonded-dict"]
        )
        # calculate dc/dh on long range
        dcdh_long = calc_dcdh(
            r,
            (g_cur_mat + g_tgt_mat) / 2,
            G_minus_g_cur_mat,
            settings["rhos"],
            settings["n_intra"],
            verbose,
        )
        dcdh_long_2D = make_matrix_2D(dcdh_long)
        # cut invert dc/dh, cut dh/dc, invert again
        dcdh_2D = cut_matrix_inverse(dcdh_long_2D, n_r, n_i, cut)
        # make it a 4D array again
        dcdh = make_matrix_4D(dcdh_2D, n_c, n_c, n_i, n_i)

    # c matrix needed for Percus-Yevick
    if settings["closure"] == "py":
        if settings["tgt_dcdh"] is not None:
            raise NotImplementedError("PY does not work with tgt_dcdh")
        k, h_hat_mat = fourier_all(r, g_cur_mat - 1)
        G_minus_g_mat = gen_interaction_matrix(
            r, input_arrays["G_minus_g_cur"], settings["non-bonded-dict"]
        )
        _, G_minus_g_hat_mat = fourier_all(r, G_minus_g_mat)
        c_mat = calc_c_matrix(
            r,
            k,
            h_hat_mat,
            G_minus_g_hat_mat,
            settings["rhos"],
            settings["n_intra"],
            verbose,
        )[cut]
    else:
        c_mat = None

    # Note: not 100% sure the order of the next three steps is correct.
    # But this is how I describe it in the paper and I tested it in for
    # neon-argon mixture and it does not make a difference

    # add the 1/g term to dc/dh and obtain inverse Jacobian
    jac_inv_mat = add_jac_inv_diagonal(
        r[cut],
        g_tgt_mat[cut],
        g_cur_mat[cut],
        dcdh[cut, cut],
        settings["rhos"],
        settings["n_intra"],
        settings["kBT"],
        settings["closure"],
        c_mat,
        verbose=verbose,
    )

    # invert the jacobian in its 2D form
    jac_mat = make_matrix_4D(
        np.linalg.inv(make_matrix_2D(jac_inv_mat)), n_c, n_c, n_i, n_i
    )
    # improve jacobian
    if settings["improve_jacobian_onset"]:
        settings["is_target_matrix"] = False
        jac_mat = improve_jacobian_onset(
            jac_mat, input_arrays, settings, verbose=verbose
        )

    # remove explicit x_ba becaue it is equal to x_ab
    jac_mat = remove_equivalent_rows_from_jacobain(jac_mat)
    return jac_mat


def remove_equivalent_rows_from_jacobain(jac_mat):
    """Averages rows and adds coluns that belong to the same interaction: x_ab, x_ba.

    Uses half-vectorization, and go row-wise through the upper triangular matrix.
    """
    n_id = jac_mat.shape[-1]
    n_t = int(n_id**0.5)
    jac_mat_reduced = jac_mat.copy()

    # add columns u_ab = u_ab + u_ba and remove u_ba
    cols_to_delete = []
    for alpha in range(n_t):
        for beta in range(n_t):
            if beta < alpha:
                index_keep = alpha * n_t + beta
                index_delete = beta * n_t + alpha
                cols_to_delete.append(index_delete)
                jac_mat_reduced[:, :, :, index_keep] += jac_mat_reduced[
                    :, :, :, index_delete
                ]
    jac_mat_reduced = np.delete(jac_mat_reduced, cols_to_delete, axis=-1)

    # average rows g_ab = (g_ab + g_ba) / 2 and remove u_ba
    rows_to_delete = []
    for alpha in range(n_t):
        for beta in range(n_t):
            if beta < alpha:
                index_keep = alpha * n_t + beta
                index_delete = beta * n_t + alpha
                rows_to_delete.append(index_delete)
                jac_mat_reduced[:, :, index_keep, :] += jac_mat_reduced[
                    :, :, index_delete, :
                ]
                jac_mat_reduced[:, :, index_keep, :] /= 2
    jac_mat_reduced = np.delete(jac_mat_reduced, rows_to_delete, axis=-2)

    return jac_mat_reduced


@if_verbose_dump_io
def calc_and_save_jacobian(input_arrays, settings, verbose=False):
    """
    Calculate Jacobian and save it

    Args:
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file
    """
    # generate it
    jac_mat = calc_jacobian(input_arrays, settings, verbose)
    # save to npz file
    if settings["out"] != "none":
        np.savez_compressed(settings["out"], jacobian=jac_mat)


@if_verbose_dump_io
def calc_dcdh(r, g_mat, G_minus_g_mat, rhos, n_intra, verbose=False):
    """
    Calculate the derivative dvec(c)/dvec(h) which is part of the Jacobian.

    Args:
        r: Distance grid
        g_mat: matrix of RDF (target or current)
        G_minus_g_cur_mat: matrix of intramolecular RDF (target or current)
        rhos: Number densities of the beads
        n_intra: Number of equal beads per molecule

    Returns:
        The inverse jacobian
    """
    # number of atom types
    n_t = len(rhos)
    # number of interactions
    n_i = int(n_t**2)
    # FT of total correlation function 'h' and G_minus_g
    k, h_hat_mat = fourier_all(r, g_mat - 1)
    _, G_minus_g_hat_mat = fourier_all(r, G_minus_g_mat)
    # Fourier matrix
    F = gen_fourier_matrix(r, fourier)
    F_inv = gen_fourier_matrix(k, fourier)
    # Ω
    Omega_hat_mat = gen_Omega_hat_mat(G_minus_g_hat_mat, rhos, n_intra)
    # ρ molecular, entries are mol densities as needed by symmetry adapted rism
    rho_mol_map = np.diag(rhos / n_intra)  # rho molecular
    # I, identity matrix
    identity = np.identity(n_t)
    # symmetry adapt h -> H
    H_hat_mat = adapt_reduced_matrix(h_hat_mat, n_intra)
    # derived from vectorizing Martin Hankes result for the derivative
    A = np.swapaxes(np.linalg.inv(Omega_hat_mat + rho_mol_map @ H_hat_mat), -1, -2)
    B = np.linalg.inv(Omega_hat_mat) @ (
        identity
        - (
            H_hat_mat
            @ np.linalg.inv(Omega_hat_mat + rho_mol_map @ H_hat_mat)
            @ rho_mol_map
        )
    )
    d_vec_C_hat_by_d_vec_H_hat = kron_2D(A, B)
    # unadapt
    d_vec_c_hat_by_d_vec_h_hat = np.zeros_like(d_vec_C_hat_by_d_vec_H_hat)
    for h, (i, j) in enumerate(itertools.product(range(n_i), range(n_i))):
        adapt_factor = np.sqrt(np.outer(n_intra, n_intra)).flatten()
        d_vec_c_hat_by_d_vec_h_hat[:, i, j] = (
            adapt_factor[j] / adapt_factor[i] * d_vec_C_hat_by_d_vec_H_hat[:, i, j]
        )
    # now it becomes an operator by diag and applying Fourier
    d_vec_c_by_d_vec_h = np.zeros((len(r), len(r), n_i, n_i))
    for h, (i, j) in enumerate(itertools.product(range(n_i), range(n_i))):
        d_vec_c_by_d_vec_h[:, :, i, j] = (
            F_inv @ np.diag(d_vec_c_hat_by_d_vec_h_hat[:, i, j]) @ F
        )
    return d_vec_c_by_d_vec_h


@if_verbose_dump_io
def calc_and_save_dcdh(input_arrays, settings, verbose=False):
    """
    Calculate dc/dh in its cut form and save it

    Args:
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file

    Returns:
        dc/dh on the cut range
    """
    # obtain r
    r = input_arrays["r"]
    # number of atom types
    n_t = len(settings["rhos"])
    # number of interactions
    n_i = int(n_t**2)
    # number of grid points in r
    n_r = len(r)
    # slices
    if settings["cut_off_res"] is None:
        cut, _ = calc_slices(r, settings["cut_off_pot"], verbose=False)
    else:
        cut, _ = calc_slices(r, settings["cut_off_res"], verbose=False)
    n_c = len(r[cut])
    # warning if short input data
    if n_c * 2 > n_r:
        print(
            "WARNING: max is smaller than twice of cut_off. This will lead to "
            "artifacts in dc/dh."
        )
    # generate matrices
    g_tgt_mat = gen_interaction_matrix(
        r, input_arrays["g_tgt"], settings["non-bonded-dict"]
    )
    # generate dc/dh, invert, cut it, and invert again
    G_minus_g_tgt_mat = gen_interaction_matrix(
        r, input_arrays["G_minus_g_tgt"], settings["non-bonded-dict"]
    )
    # calculate dc/dh on long range
    dcdh_long = calc_dcdh(
        r,
        g_tgt_mat,
        G_minus_g_tgt_mat,
        settings["rhos"],
        settings["n_intra"],
        verbose,
    )
    dcdh_long_2D = make_matrix_2D(dcdh_long)
    # cut invert dc/dh, cut dh/dc, invert again
    dcdh_2D = cut_matrix_inverse(dcdh_long_2D, n_r, n_i, cut)
    # make it a 4D array again
    dcdh = make_matrix_4D(dcdh_2D, n_c, n_c, n_i, n_i)
    # save to npz file
    if settings["out"] != "none":
        np.savez_compressed(settings["out"], dcdh=dcdh)


@if_verbose_dump_io
def add_jac_inv_diagonal(
    r,
    g_tgt_mat,
    g_cur_mat,
    dcdh,
    rhos,
    n_intra,
    kBT,
    closure,
    c_mat,
    verbose=False,
):
    """
    Calculate du/dg, the inverse of the Jacobian.

    Args:
        r: Distance grid
        g_tgt_mat: target RDFs
        g_cur_mat: current RDFs
        dcdh_2D: derivative dc/dh
        rhos: Number densities of the beads
        n_intra: Number of equal beads per molecule
        kBT: Boltzmann constant times temperature
        closure: OZ-equation closure ('hnc' or 'py')

    Returns:
        Matrix inverse of the Jacobian
    """
    # number of atom types
    n_t = len(rhos)
    # number of interactions
    n_i = int(n_t**2)
    # vectorize RDF matrices
    g_tgt_vec = vectorize_full(g_tgt_mat)
    g_cur_vec = vectorize_full(g_cur_mat)
    # average RDF for better stability
    g_avg_vec = (g_tgt_vec + g_cur_vec) / 2
    jac_inv_mat = np.zeros((len(r), len(r), n_i, n_i))
    if closure == "hnc":
        with np.errstate(divide="ignore", invalid="ignore", under="ignore"):
            for i, j in itertools.product(range(n_i), range(n_i)):
                if i == j:
                    diagonal_term = np.diag(1 - 1 / g_avg_vec[:, i])
                else:
                    diagonal_term = 0
                jac_inv_mat[:, :, i, j] = kBT * (diagonal_term - dcdh[:, :, i, j])
    elif closure == "py":
        # implementation is a bit tricky because one needs c_mat again, but
        # if tgt_dcdh the intramolecular interactions are not present
        c_vec = vectorize_full(c_mat)
        with np.errstate(divide="ignore", invalid="ignore", under="ignore"):
            for i, j in itertools.product(range(n_i), range(n_i)):
                if i == j:
                    same_nb_term = np.diag(c_vec[:, i] / g_avg_vec[:, i])
                else:
                    same_nb_term = 0
                jac_inv_mat[:, :, i, j] = kBT * (
                    (same_nb_term - dcdh[:, :, i, j]) / (g_avg_vec[:, i] - c_vec[:, i])
                )
                # times g version (does not work)
                # jac_inv_mat[:, :, i, j] = kBT * (
                # (same_nb_term - g_avg_vec[:, i] * dcdh[:, :, i, j])
                # / (g_avg_vec[:, i] ** 2 - g_avg_vec[:, i] * c_vec[:, i])
                # )
    # make negative infinities a large negative number, increasing stability
    jac_inv_mat[np.isneginf(jac_inv_mat)] = -1e37
    return jac_inv_mat


def gen_Omega_hat_mat(G_minus_g_hat_mat, rhos, n_intra):
    # σ is any row sum of ω
    sigma_R = G_minus_g_hat_mat @ np.diag(rhos) + np.identity(len(rhos))
    # weighting of row sum σ
    Omega_hat_mat = np.diag(np.sqrt(n_intra)) @ sigma_R @ np.diag(1 / np.sqrt(n_intra))
    return Omega_hat_mat


def adapt_reduced_matrix(mat, n_intra):
    """Adapt the prefactors of a matrix to be compatible with the symmetry reduced RISM
    equation.

    The input matrix is already reduced (rows are atom types not atoms), but factors
    need to be applied."""
    Mat = np.diag(np.sqrt(n_intra)) @ mat @ np.diag(np.sqrt(n_intra))
    return Mat


def unadapt_reduced_matrix(Mat, n_intra):
    """Unadapt the prefactors of a matrix compatible with the symmetry reduced RISM
    equation back to the regular form.

    The input matrix is already reduced (rows are atom types not atoms) and adapted.
    Factors are applied to get back to the regular matrix."""
    mat = np.diag(1 / np.sqrt(n_intra)) @ Mat @ np.diag(1 / np.sqrt(n_intra))
    return mat


if __name__ == "__main__":
    main()
