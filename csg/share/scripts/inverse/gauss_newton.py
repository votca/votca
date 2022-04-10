#!/usr/bin/env python3
"""Multi purpose script for updating Potentials with Newton or Gauss-Newton method."""
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
# Delta_U: potential update (u_{k+1} - u_k) = Δu
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
import numbers
import numpy as np
import sys
import xml.etree.ElementTree as ET

from csg_functions import (
    calc_grid_spacing,
    calc_slices,
    devectorize,
    eval_expr,
    extrapolate_Delta_u_left_constant,
    gen_beadtype_property_array,
    gen_flag_isfinite,
    gen_interaction_dict,
    gen_interaction_matrix,
    get_bead_types,
    get_density_dict,
    get_n_intra_dict,
    get_non_bonded,
    if_verbose_dump_io,
    kron_2D,
    make_matrix_2D,
    read_all_tables,
    save_tables,
    solve_lsq_with_linear_constraints,
    triangular_number,
    vectorize,
)

if not sys.version_info >= (3, 5):
    raise Exception("This script needs Python 3.6+.")

# constants
BAR_PER_MD_PRESSURE = 16.6053904  # md pressure is kJ/mol/nm³ as in Gromacs

# raise all numpy errors. If errors are expected, use np.errstate()
np.seterr(all="raise")


def main():
    # get command line arguments
    args = get_args()
    # process and prepare input
    input_arrays, settings = process_input(args)
    # guess potential from distribution
    # newton update
    if settings["subcommand"] == "newton":
        output_arrays = newton_update(
            input_arrays, settings, verbose=settings["verbose"]
        )
    # gauss-newton update
    elif settings["subcommand"] == "gauss-newton":
        if settings["multistate"]:
            output_arrays = multistate_gauss_newton_update(
                input_arrays, settings, verbose=settings["verbose"]
            )
        else:
            if any(settings["upd_pots"].values()):  # only update if requested
                output_arrays = gauss_newton_update(
                    input_arrays, settings, verbose=settings["verbose"]
                )
            else:
                print("No potentials to update with iie gauss-newton, writing zeros.")
                output_arrays = zero_update(
                    input_arrays, settings, verbose=settings["verbose"]
                )
    # save all potential updates to table files
    save_tables(output_arrays, settings)


def get_args(iie_args=None):
    """Define and parse command line arguments.

    If iie_args is given, parse them instead of cmdlineargs.
    """
    description = "Calculate u or Δu with Integral Equations."
    parser = argparse.ArgumentParser(description=description)
    # subparsers
    subparsers = parser.add_subparsers(dest="subcommand")
    parser_newton = subparsers.add_parser(
        "newton", help="potential update using Newton method"
    )
    parser_gauss_newton = subparsers.add_parser(
        "gauss-newton", help="potential update using Gauss-Newton method"
    )
    # all subparsers
    for pars in (parser_newton, parser_gauss_newton):
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
            "--jacobian",
            type=argparse.FileType("r"),
            nargs="+",
            default=None,
            help=".npz file with the Jacobian. Multiple files for multistate.",
        ),
        pars.add_argument(
            "--volume",
            type=float,
            required=True,
            nargs="+",
            metavar="VOL",
            help="the volume of the box. Multiple if multistate",
        )
        pars.add_argument(
            "--topol",
            type=argparse.FileType("r"),
            required=True,
            nargs="+",
            metavar="TOPOL",
            help="XML topology file. Multiple if multistate",
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
            "--g-cur-ext",
            type=str,
            required=True,
            metavar="RDF_CUR_EXT",
            help="extension of current RDF files",
        )
        pars.add_argument(
            "--flatten-at-cut-off",
            action="store_const",
            const=True,
            default=False,
            help=(
                "Weather to set the last point of Δu to zero. "
                "The pressure constraint adapts to one point less."
            ),
        )
    # GN only options
    parser_gauss_newton.add_argument(
        "--pressure-constraint",
        type=str,
        default=None,
        dest="pressure_constraint",
        nargs="*",
        help=(
            'String of form ",p_tgt,p_cur". Starting '
            "comma is needed to prevent confusion when "
            "p_tgt is negative"
        ),
    )
    parser_gauss_newton.add_argument(
        "--kirkwood-buff-constraint",
        action="store_const",
        const=True,
        default=False,
        help="Whether to match the KBI of the target RDF",
    )
    parser_gauss_newton.add_argument(
        "--residual-weighting",
        dest="residual_weighting",
        type=str,
        required=True,
    )
    parser_gauss_newton.add_argument(
        "--tgt-dists",
        type=str,
        required=True,
        metavar="TGT_DISTS",
        help=(
            "Which distributions are targeted. Pairs of names and "
            "bool. Pair separated by comma. Pairs separated by colon. "
            "true for target, false for not target."
        ),
    )
    parser_gauss_newton.add_argument(
        "--upd-pots",
        type=str,
        required=True,
        metavar="UPD_POTS",
        help=(
            "Which potentials are to be modified. Pairs of names and "
            "numbers. Pair separated by comma. Pairs separated by colon. "
            "1 for update, 0 for no update."
        ),
    )
    # GN and dcdh subparsers
    parser_gauss_newton.add_argument(
        "--multistate",
        dest="multistate",
        action="store_const",
        const=True,
        default=False,
        help="enable multistate method",
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
    # multistate settings
    multistate = False
    if args.subcommand == "gauss-newton":
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
        "flatten_at_cut_off",
        "out",
        "residual_weighting",
        "subcommand",
        "verbose",
    )
    settings = {key: vars(args)[key] for key in args_to_copy if key in vars(args)}
    settings["non-bonded-dict"] = non_bonded_dict
    settings["rhos"] = rhos
    settings["n_intra"] = n_intra
    settings["r0-removed"] = r0_removed

    # load Joacbian
    settings["jacobians"] = []  # one per state
    for jac_file in args.jacobian:
        try:
            # close file(s), because we use np.load on file name
            jac_file.close()
            settings["jacobians"].append(np.load(jac_file.name)["jacobian"])
        except (FileNotFoundError, ValueError):
            raise Exception("Can not load jacobian file(s) that were provided")
    if multistate:
        assert len(settings["jacobians"]) == len(state_names)
    else:
        assert len(settings["jacobians"]) == 1

    # determine cut-off xml path
    cut_off_path = {
        "newton": "./inverse/newton/cut_off",
        "gauss-newton": "./inverse/gauss_newton/cut_off",
    }[args.subcommand]
    # try to read it from xml ET
    try:
        settings["cut_off_pot"] = float(options.find(cut_off_path).text)
    except (AttributeError, ValueError):
        raise Exception(cut_off_path + " must be a float in settings.xml")
    # determine cut-off for residuals
    if args.subcommand == "gauss-newton":
        cut_residual = options.find("./inverse/gauss_newton/cut_residual")
        settings["cut_off_res"] = float(cut_residual.text)

    # constraints
    if args.subcommand == "gauss-newton":
        constraints = []
        if args.pressure_constraint is not None:
            if multistate:
                raise NotImplementedError("multistate + constraints not implemented")
            else:
                p_target = float(args.pressure_constraint[0].lstrip(",").split(",")[0])
                p_current = float(args.pressure_constraint[0].lstrip(",").split(",")[1])
                constraints.append(
                    {"type": "pressure", "target": p_target, "current": p_current}
                )
        if args.kirkwood_buff_constraint:
            constraints.append({"type": "kirkwood-buff-integral"})
        settings["constraints"] = constraints

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

    # which distributions to target and which potentials to update
    if args.subcommand == "gauss-newton":
        settings["tgt_dists"] = {
            pair[0]: pair[1].lower().strip() == "true"
            for pair in [
                pair.split(",") for pair in args.tgt_dists.strip(":").split(":")
            ]
        }
        settings["upd_pots"] = {
            pair[0]: pair[1].strip() == "1"
            for pair in [
                pair.split(",") for pair in args.upd_pots.strip(":").split(":")
            ]
        }
        # check that there is a value for each non-bonded interaction
        assert set(settings["tgt_dists"].keys()) == set(non_bonded_dict.keys())
        assert set(settings["upd_pots"].keys()) == set(non_bonded_dict.keys())
        # check that not all tgt_dist are False
        if not any(settings["tgt_dists"].values()):
            raise ValueError("No distribution set to be target.")

    return input_arrays, settings


@if_verbose_dump_io
def newton_update(input_arrays, settings, verbose=False):
    """Calculate Newton potential update.

    Args:
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file

    Returns:
        dictionary of potential updates including flags to be saved
    """
    # obtain r
    r = input_arrays["r"]
    # number of atom types
    n_t = len(settings["rhos"])
    # number of interactions including redundand ones
    n_i = triangular_number(n_t)
    # slices
    cut, tail = calc_slices(r, settings["cut_off_pot"], verbose=False, offset=-1)
    n_c = len(r[cut])
    # generate matrices and directly cut them
    g_tgt_mat = gen_interaction_matrix(
        r, input_arrays["g_tgt"], settings["non-bonded-dict"]
    )[cut]
    g_cur_mat = gen_interaction_matrix(
        r, input_arrays["g_cur"], settings["non-bonded-dict"]
    )[cut]
    # get Jacobian and cut away flesh
    # (usually cuts just one row/col, since we do not update last point of potential)
    jac_mat = settings["jacobians"][0][cut, cut, :, :]
    # make it 2D
    jac_mat_2D = make_matrix_2D(jac_mat)
    # Delta g for potential update
    Delta_g_mat = g_cur_mat - g_tgt_mat
    # vectorize Delta g
    Delta_g_vec = vectorize(Delta_g_mat)
    # flatten Delta g
    Delta_g_flat = Delta_g_vec.T.flatten()

    # Newton update in the flat/2D representation
    with np.errstate(invalid="ignore"):
        Delta_u_flat = -1 * np.linalg.lstsq(jac_mat_2D, Delta_g_flat, rcond=None)[0]

    # make interactions an own dimension again
    Delta_u_vec = Delta_u_flat.reshape((n_i, n_c)).T
    # Δu matrix
    Delta_u_mat = devectorize(Delta_u_vec)
    # prepare output
    output_arrays = {}
    for non_bonded_name, Delta_u_dict in gen_interaction_dict(
        r[cut], Delta_u_mat, settings["non-bonded-dict"]
    ).items():
        Delta_u = Delta_u_dict["y"]
        Delta_u_flag = gen_flag_isfinite(Delta_u)
        # add value at r=0 if it was removed
        if settings["r0-removed"]:
            r_out = np.concatenate(([0.0], r[cut]))
            Delta_u = np.concatenate(([np.nan], Delta_u))
            Delta_u_flag = np.concatenate((["o"], Delta_u_flag))
        else:
            r_out = r[cut]
        # change NaN in the core region to first valid value
        Delta_u = extrapolate_Delta_u_left_constant(Delta_u, Delta_u_flag)
        # save for output
        output_arrays[non_bonded_name] = {
            "x": r_out,
            "y": Delta_u,
            "flag": Delta_u_flag,
        }
    return output_arrays


def gen_indices_upd_pots_tgt_dists(settings):
    index_tgt_dists = []
    index_not_tgt_dists = []
    index_upd_pots = []
    index_not_upd_pots = []
    bead_types = get_bead_types(settings["non-bonded-dict"])
    non_bonded_dict_inv = {v: k for k, v in settings["non-bonded-dict"].items()}
    for i, ((b1, bead1), (b2, bead2)) in enumerate(
        (
            ((b1, bead1), (b2, bead2))
            for b1, bead1 in enumerate(bead_types)
            # due to row-wise half-vectorization
            for b2, bead2 in [
                (b2, bead2) for (b2, bead2) in enumerate(bead_types) if b2 >= b1
            ]
        )
    ):
        interaction_name = non_bonded_dict_inv[frozenset({bead1, bead2})]
        if settings["tgt_dists"][interaction_name] is True:
            index_tgt_dists.append(i)
        else:
            index_not_tgt_dists.append(i)
        if settings["upd_pots"][interaction_name] is True:
            index_upd_pots.append(i)
        else:
            index_not_upd_pots.append(i)
    return (
        index_tgt_dists,
        index_not_tgt_dists,
        index_upd_pots,
        index_not_upd_pots,
    )


def gen_dfdu_matrix(n_c_pot, n_c_res, Delta_r):
    """Generate with df/du partial derivatives.

    There are as many grid points for forces as there are for the potential (w/o
    cut-off). The grid of the forces is shifted relative to the potential/RDF grid:
    r + Δr/2.
    """
    dfdu = np.zeros((n_c_pot, n_c_res))
    # create indices for changing the main- and super-diagonal
    main_diag = np.diag_indices(n_c_pot, ndim=2)
    # copy needed, because diag_indices generates twice the same object
    super_diag = list(np.diag_indices(n_c_pot - 1, ndim=2)[i].copy() for i in (0, 1))
    # shift to geth the super diagonal
    super_diag[1] += 1
    super_diag = tuple(super_diag)  # tuple needed for indexing
    dfdu[main_diag] = -1 / Delta_r
    dfdu[super_diag] = +1 / Delta_r
    return dfdu


@if_verbose_dump_io
def gauss_newton_update(input_arrays, settings, verbose=False):
    """Calculate Gauss-Newton potential update based on s.a. RISM-OZ and closure.

    Args:
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file

    Returns:
        dictionary of potential updates including flags to be saved
    """
    # obtain r
    r = input_arrays["r"]
    # number of atom types
    n_t = len(settings["rhos"])
    # number of interactions including redundand ones
    n_i = triangular_number(n_t)
    # grid spacing
    Delta_r = calc_grid_spacing(r)
    # slices and indices
    # we don't want to update at cut-off, therefore offset=-1
    cut_pot, tail_pot = calc_slices(
        r, settings["cut_off_pot"], verbose=False, offset=-1
    )
    # sometimes we need to refer to values including the value of the cut-off
    cut_pot_p1, tail_pot_p1 = calc_slices(
        r, settings["cut_off_pot"], verbose=False, offset=0
    )
    # residuum cut also -= 1 for consistency when cut-residuum = cut-off
    cut_res, tail_res = calc_slices(
        r, settings["cut_off_res"], verbose=False, offset=-1
    )
    n_c_pot = len(r[cut_pot])
    n_c_res = len(r[cut_res])

    # generate matrices from dicts
    g_tgt_mat = gen_interaction_matrix(
        r, input_arrays["g_tgt"], settings["non-bonded-dict"]
    )
    g_cur_mat = gen_interaction_matrix(
        r, input_arrays["g_cur"], settings["non-bonded-dict"]
    )
    # get Jacobian
    jac_mat = settings["jacobians"][0][cut_res, cut_res, :, :]
    # weighting
    g_tgt_vec = vectorize(g_tgt_mat)
    g_cur_vec = vectorize(g_cur_mat)
    weights = gen_residual_weights(
        settings["residual_weighting"],
        r[cut_res],
        n_c_res,
        n_i,
        1,
        g_tgt_vec[cut_res],
        [1],
    )
    # Delta g for potential update
    Delta_g_mat = g_cur_mat - g_tgt_mat
    # vectorize Delta g
    Delta_g_vec = vectorize(Delta_g_mat)

    # one might not want to take all dists in account or update all potentials
    # determine which rows/cols should be removed
    (
        index_tgt_dists,
        index_not_tgt_dists,
        index_upd_pots,
        index_not_upd_pots,
    ) = gen_indices_upd_pots_tgt_dists(settings)
    n_tgt_dists = len(index_tgt_dists)
    n_upd_pots = len(index_upd_pots)
    # remove rows and columns from jacobian
    jac_mat = np.delete(jac_mat, index_not_tgt_dists, axis=-2)  # del rows
    jac_mat = np.delete(jac_mat, index_not_upd_pots, axis=-1)  # del cols
    # remove cols from weights
    weights = np.delete(weights, index_not_tgt_dists, axis=1)  # del cols
    # remove cols from Δg vector
    Delta_g_vec = np.delete(Delta_g_vec, index_not_tgt_dists, axis=1)  # del cols
    # remove rows and columns from target RDF matrix for pressure constraint
    g_tgt_vec = vectorize(g_tgt_mat)
    g_tgt_vec = np.delete(g_tgt_vec, index_not_upd_pots, axis=-2)  # del rows
    g_tgt_vec = np.delete(g_tgt_vec, index_not_upd_pots, axis=-1)  # del cols

    # make Jacobian 2D
    jac_2D = make_matrix_2D(jac_mat)
    # prepare A0 and C matrix
    # update potential, A0 just cuts of the tail
    A0_2D = kron_2D(
        np.identity(n_upd_pots),
        np.append(np.identity(n_c_pot), np.zeros((n_c_res - n_c_pot, n_c_pot)), axis=0),
    )

    # weight Jacobian
    jac_2D_weighted = np.diag(weights.T.flatten()) @ jac_2D
    # cut and weight Delta_g and obtain residuals
    residuals_vec = np.zeros((n_c_res, n_tgt_dists))
    for i in range(n_tgt_dists):
        residuals_vec[:, i] = np.diag(weights[:, i]) @ Delta_g_vec[cut_res, i]
    # flatten residuals
    residuals_flat = residuals_vec.T.flatten()
    # we will solve the least squares problem |A @ x - b| with constraints C @ x - d
    A = jac_2D_weighted @ A0_2D
    b = residuals_flat
    # nr of constraints
    n_constraints_rows = {
        "pressure": 1,
        "kirkwood-buff-integral": n_tgt_dists,
        "potential-energy": 1,
    }
    n_constraints = sum(
        n_constraints_rows[constr["type"]] for constr in settings["constraints"]
    )
    d = np.zeros(n_constraints)
    C = np.zeros((n_constraints, n_upd_pots * n_c_pot))
    # prepare all constraints
    c_run = 0  #
    for c, constraint in enumerate(settings["constraints"]):

        if constraint["type"] == "pressure":
            # current pressure
            p = constraint["current"] / BAR_PER_MD_PRESSURE
            # target pressure
            p_tgt = constraint["target"] / BAR_PER_MD_PRESSURE
            # print if verbose
            if settings["verbose"]:
                print(
                    f"Constraining pressure: target is {constraint['target']} bar, "
                    f"current value is {constraint['current']} bar"
                )
            # Written in matrix form but equivalent to paper
            # define df/du matrix
            # use short range for both f and u
            dfdu = gen_dfdu_matrix(n_c_pot, n_c_pot, Delta_r)
            dfdu_all = np.kron(np.eye(n_upd_pots), dfdu)
            # volume of sphere shell element
            r3_dr = np.tile(
                ((r[cut_pot] + Delta_r) ** 4 - r[cut_pot] ** 4) / 4, n_upd_pots
            )
            # average RDF in the shell
            g_cur_avg = (
                g_cur_vec[cut_pot_p1][:-1].T.flatten()
                + g_cur_vec[cut_pot_p1][1:].T.flatten()
            ) / 2
            # density product ρ_i * ρ_j as vector of same length as r_i
            rho_ab = np.repeat(
                vectorize(np.outer(*([settings["rhos"]] * 2)))[index_upd_pots], n_c_pot
            )
            # extra factor of 2 for the off-diagonal interactions
            # this factor arises from the double sum over bead types α and β
            extra_mat = np.ones((n_t, n_t)) * 2
            np.fill_diagonal(extra_mat, 1)
            extra_factor = np.repeat(vectorize(extra_mat)[index_upd_pots], n_c_pot)
            # dp/df
            dpdf = -2 / 3 * np.pi * rho_ab * g_cur_avg * r3_dr * extra_factor
            C[c, :] = dpdf @ dfdu_all
            if settings["flatten_at_cut_off"]:
                C[
                    c, n_c_pot - 1 :: n_c_pot
                ] = 0  # no constraint for last point of each Δu
            d[c] = p - p_tgt

        elif constraint["type"] == "kirkwood-buff-integral":
            for i in range(n_tgt_dists):
                # we leave out pre factors
                # current KBI
                G = (
                    4
                    * np.pi
                    * np.sum(r[cut_res] ** 2 * (g_cur_vec[cut_res, i] - 1))
                    * Delta_r
                )
                # target KBI
                # TODO: custom target not implemented
                # if "target" in constraint:
                # G_tgt = constraint["target"]
                G_tgt = (
                    4
                    * np.pi
                    * np.sum(r[cut_res] ** 2 * (g_tgt_vec[cut_res, i] - 1))
                    * Delta_r
                )
                # print if verbose
                if settings["verbose"]:
                    print(
                        f"Constraining Kirkwood-Buff integral of interaction {i} from "
                        f"0 to r={r[cut_res][-1]:.4f} nm: target is {G_tgt:.5f} nm³, "
                        f"current value is {G:.5f} nm³"
                    )
                # define C row
                dGdg = 4 * np.pi * Delta_r * r[cut_res] ** 2
                C[c_run + i, :] = dGdg @ make_matrix_2D(
                    jac_mat[cut_res, cut_pot, i : i + 1, :]
                )
                if settings["flatten_at_cut_off"]:
                    C[
                        c_run + i, n_c_pot - 1 :: n_c_pot
                    ] = 0  # no constraint for last point of each Δu
                d[c_run + i] = G - G_tgt

        elif constraint["type"] == "potential-energy":
            # not yet in use
            if n_t > 1:
                raise NotImplementedError(
                    "potential energy constraint not implemented "
                    "for more than one bead"
                )
            # we leave out pre factors (rho should come back for multiple beads)
            # current PE
            PE = constraint["current"]
            # target PE
            PE_tgt = constraint["target"]
            # print if verbose
            if settings["verbose"]:
                print(
                    f"Constraining potential energy: target is "
                    f"{constraint['target']} kJ/mol, "
                    f"current value is {constraint['current']} kJ/mol"
                )
            rho_i = np.repeat(
                vectorize(np.outer(*([settings["rhos"]] * 2)))[index_upd_pots], n_c_pot
            )
            # define C row
            C[c, :] = (
                2 * np.pi * rho_i * r[cut_pot] ** 2 * g_tgt_vec[cut_pot].T.flatten()
            )
            if settings["flatten_at_cut_off"]:
                C[
                    c, n_c_pot - 1 :: n_c_pot
                ] = 0  # no constraint for last point of each Δu
            d[c] = PE - PE_tgt
        else:
            raise NotImplementedError(
                "not implemented constraint type: " + constraint["type"]
            )
        # increal running index
        c_run += n_constraints_rows[constraint["type"]]
    # solve least squares problem with linear constraints
    # if C.shape[0] == 0 and A.shape[0] == A.shape[1]:
    # dx_flat = np.linalg.solve(A, b)
    dx_flat = solve_lsq_with_linear_constraints(A, C, b, d)
    # obtain potential update
    Delta_u_flat = -A0_2D @ dx_flat
    Delta_u_vec = Delta_u_flat.reshape(n_upd_pots, n_c_res).T
    # insert zeros for ignored potentials
    for i in index_not_upd_pots:
        Delta_u_vec = np.insert(Delta_u_vec, i, np.zeros(n_c_res), axis=1)
    # Delta_u matrix
    Delta_u_mat = devectorize(Delta_u_vec)
    # cut off tail
    # Delta_u_mat = Delta_u_mat[cut_pot]
    # prepare output
    output_arrays = {}
    for non_bonded_name, Delta_u_dict in gen_interaction_dict(
        r[cut_res], Delta_u_mat, settings["non-bonded-dict"]
    ).items():
        r_out = Delta_u_dict["x"]
        Delta_u = Delta_u_dict["y"]
        Delta_u_flag = gen_flag_isfinite(Delta_u)
        Delta_u_flag[tail_pot_p1] = "o"
        # shift potential to make value before cut-off zero
        if settings["flatten_at_cut_off"]:
            # hmmmmmm will probably throw this out
            Delta_u[cut_pot] -= Delta_u[cut_pot][-1]
            # Delta_u[cut_pot][-1] = 0  # this one does not really make it flat
        # add value at r=0 if it was removed
        if settings["r0-removed"]:
            r_out = np.concatenate(([0.0], r_out))
            Delta_u = np.concatenate(([np.nan], Delta_u))
            Delta_u_flag = np.concatenate((["o"], Delta_u_flag))
        # change NaN in the core region to first valid value
        Delta_u = extrapolate_Delta_u_left_constant(Delta_u, Delta_u_flag)
        # save for output
        output_arrays[non_bonded_name] = {
            "x": r_out,
            "y": Delta_u,
            "flag": Delta_u_flag,
        }
    return output_arrays


@if_verbose_dump_io
def multistate_gauss_newton_update(input_arrays, settings, verbose=False):
    """Calculate Gauss-Newton potential update based on s.a. RISM-OZ and closure.

    Args:
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file

    Returns:
        dictionary of potential updates including flags to be saved
    """
    # ######## MULTISTATE ######### #
    # MULTISTATE ################## #
    # ################## MULTISTATE #
    # ######## MULTISTATE ######### #

    # obtain r
    r = input_arrays["r"]
    # number of atom types
    n_t = len(settings["rhos"][0])
    # number of interactions including redundand ones
    n_i = triangular_number(n_t)
    # number of states
    state_names = settings["state_names"]  # shortcut
    n_s = len(state_names)
    # slices
    # we don't want to update at cut-off, therefore offset=-1
    cut_pot, tail_pot = calc_slices(
        r, settings["cut_off_pot"], verbose=False, offset=-1
    )
    # sometimes we need to refer to values including the value of the cut-off
    cut_pot_p1, tail_pot_p1 = calc_slices(
        r, settings["cut_off_pot"], verbose=False, offset=0
    )
    cut_res, tail_res = calc_slices(
        r, settings["cut_off_res"], verbose=False, offset=-1
    )
    n_c_pot = len(r[cut_pot])
    n_c_res = len(r[cut_res])
    # generate matrices
    # g_tgt and g_cur, [r, state * interaction]
    g_tgt_vec = np.concatenate(
        [
            vectorize(
                gen_interaction_matrix(
                    r, input_arrays[s]["g_tgt"], settings["non-bonded-dict"]
                )
            )
            for s in state_names
        ],
        axis=1,
    )
    g_cur_vec = np.concatenate(
        [
            vectorize(
                gen_interaction_matrix(
                    r, input_arrays[s]["g_cur"], settings["non-bonded-dict"]
                )
            )
            for s in state_names
        ],
        axis=1,
    )
    # create the jacobian
    jac_mat = np.concatenate(settings["jacobians"], axis=2)[cut_res, cut_res, :, :]
    # make Jacobian 2D
    jac_2D = make_matrix_2D(jac_mat)
    # weights
    weights = gen_residual_weights(
        settings["residual_weighting"],
        r[cut_res],
        n_c_res,
        n_i,
        n_s,
        g_tgt_vec[cut_res],
        settings["state_weights"],
    )
    # weight Jacobian
    jac_2D = np.diag(weights.T.flatten()) @ jac_2D
    # Delta g for potential update
    Delta_g_vec = g_cur_vec - g_tgt_vec
    # cut and weight Delta_g and obtain residuals
    residuals_vec = np.zeros((n_c_res, n_s * n_i))
    for s in range(n_s):
        for i in range(n_i):
            residuals_vec[:, s * n_i + i] = (
                np.diag(weights[:, s * n_i + i]) @ Delta_g_vec[cut_res, s * n_i + i]
            )
    # flatten residuals
    residuals_flat = residuals_vec.T.flatten()
    # update potential, A0 just cuts of the tail
    A0_2D = kron_2D(
        np.identity(n_i),
        np.append(np.identity(n_c_pot), np.zeros((n_c_res - n_c_pot, n_c_pot)), axis=0),
    )
    A = jac_2D @ A0_2D
    b = residuals_flat
    # no constraints
    # C = np.zeros((0, n_i * n_c_pot))
    # d = np.zeros(0)
    # solve linear equation
    # Delta_u_flat = -A0_2D @ solve_lsq_with_linear_constraints(A, C, b, d)
    Delta_u_flat = -A0_2D @ np.linalg.lstsq(A, b, rcond=None)[0]
    Delta_u_vec = Delta_u_flat.reshape(n_i, n_c_res).T
    # Delta_u matrix
    Delta_u_mat = devectorize(Delta_u_vec)
    # Delta_u_mat = Delta_u_mat[cut_pot]
    # prepare output
    output_arrays = {}
    for non_bonded_name, Delta_u_dict in gen_interaction_dict(
        r[cut_res], Delta_u_mat, settings["non-bonded-dict"]
    ).items():
        r_out = Delta_u_dict["x"]
        Delta_u = Delta_u_dict["y"]
        Delta_u_flag = gen_flag_isfinite(Delta_u)
        Delta_u_flag[tail_pot_p1] = "o"
        # shift potential to make value before cut-off zero
        if settings["flatten_at_cut_off"]:
            Delta_u[cut_pot] -= Delta_u[cut_pot][-1]
        # add value at r=0 if it was removed
        if settings["r0-removed"]:
            r_out = np.concatenate(([0.0], r_out))
            Delta_u = np.concatenate(([np.nan], Delta_u))
            Delta_u_flag = np.concatenate((["o"], Delta_u_flag))
        else:
            r_out = r
        # change NaN in the core region to first valid value
        Delta_u = extrapolate_Delta_u_left_constant(Delta_u, Delta_u_flag)
        # save for output
        output_arrays[non_bonded_name] = {
            "x": r_out,
            "y": Delta_u,
            "flag": Delta_u_flag,
        }
    return output_arrays


@if_verbose_dump_io
def zero_update(input_arrays, settings, verbose=False):
    """Calculate Gauss-Newton potential update based on s.a. RISM-OZ and closure.

    Args:
        input_arrays: nested dict holding the distributions
        settings: dict holding relevant settings
        verbose: save parameters and return of this and contained functions as numpy
                 file

    Returns:
        dictionary of potential updates including flags to be saved
    """
    # obtain r
    r = input_arrays["r"]
    # number of atom types
    n_t = len(settings["rhos"])
    # slices
    # we don't want to update at cut-off, therefore offset=-1
    cut_pot, tail_pot = calc_slices(
        r, settings["cut_off_pot"], verbose=False, offset=-1
    )
    # sometimes we need to refer to values including the value of the cut-off
    _, tail_pot_p1 = calc_slices(r, settings["cut_off_pot"], verbose=False, offset=0)
    cut_res, tail_res = calc_slices(r, settings["cut_off_res"], verbose=False)
    n_c_res = len(r[cut_res])
    Delta_u_mat = np.zeros((n_c_res, n_t, n_t))
    # prepare output
    output_arrays = {}
    for non_bonded_name, Delta_u_dict in gen_interaction_dict(
        r[cut_res], Delta_u_mat, settings["non-bonded-dict"]
    ).items():
        r_out = Delta_u_dict["x"]
        Delta_u = Delta_u_dict["y"]
        Delta_u_flag = gen_flag_isfinite(Delta_u)
        Delta_u_flag[tail_pot_p1] = "o"
        # add value at r=0 if it was removed
        if settings["r0-removed"]:
            r_out = np.concatenate(([0.0], r_out))
            Delta_u = np.concatenate(([0.0], Delta_u))
            Delta_u_flag = np.concatenate((["o"], Delta_u_flag))
        # save for output
        output_arrays[non_bonded_name] = {
            "x": r_out,
            "y": Delta_u,
            "flag": Delta_u_flag,
        }
    return output_arrays


def gen_residual_weights(
    res_weight_scheme, r, n_c_res, n_i, n_s, g_tgt_vec, state_weights
):
    """Generate weights with one column along r for each interaction."""
    with np.errstate(all="ignore"):
        weights = eval_expr(
            res_weight_scheme.strip("'"),
            variables={
                "r": np.repeat(np.transpose([r]), n_s * n_i, 1),
                "g_tgt": g_tgt_vec,
                "max_r": max(r),
            },
        )
    # if weights is a float, all weights will be the same
    if isinstance(weights, numbers.Number):
        weights = np.ones((n_c_res, n_s * n_i)) * weights
    elif weights.shape == (n_c_res, n_s * n_i):
        pass
    else:
        raise Exception(
            "Something went wrong: weights should be float or a np.array "
            f"of shape ({n_c_res}, {n_s * n_i}). Instead it is {weights}"
        )

    # weight each state by state_weigths
    for s in range(n_s):
        weights[:, n_i * s : n_i * (s + 1)] *= state_weights[s]

    return weights


if __name__ == "__main__":
    main()
