/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/bse.h"
#include "votca/xtp/ecpbasisset.h"
#include "votca/xtp/gwbse.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/vxc_potential.h"

using boost::format;
using namespace boost::filesystem;
using std::flush;
namespace votca {
namespace xtp {

Index GWBSE::CountCoreLevels() {
  Index ignored_corelevels = 0;
  if (!_orbitals.hasECPName()) {
    ECPBasisSet basis;
    basis.Load("corelevels");
    Index coreElectrons = 0;
    for (const auto& atom : _orbitals.QMAtoms()) {
      coreElectrons += basis.getElement(atom.getElement()).getNcore();
    }
    ignored_corelevels = coreElectrons / 2;
  }
  return ignored_corelevels;
}

void GWBSE::Initialize(tools::Property& options) {

  std::string key = Identify();

  // getting level ranges
  Index rpamax = 0;
  Index rpamin = 0;  // never changes
  Index qpmin = 0;
  Index qpmax = 0;
  Index bse_vmin = 0;
  Index bse_cmax = 0;

  Index homo = _orbitals.getHomo();  // indexed from 0
  Index num_of_levels = _orbitals.getBasisSetSize();
  Index num_of_occlevels = _orbitals.getNumberOfAlphaElectrons();

  std::string ranges = options.get(key + ".ranges").as<std::string>();

  // now check validity, and get rpa, qp, and bse level ranges accordingly

  if (ranges == "factor") {

    double rpamaxfactor = options.get(key + ".rpamax").as<double>();
    rpamax = Index(rpamaxfactor * double(num_of_levels)) -
             1;  // total number of levels

    double qpminfactor = options.get(key + ".qpmin").as<double>();
    qpmin =
        num_of_occlevels - Index(qpminfactor * double(num_of_occlevels)) - 1;

    double qpmaxfactor = options.get(key + ".qpmax").as<double>();
    qpmax =
        num_of_occlevels + Index(qpmaxfactor * double(num_of_occlevels)) - 1;

    double bseminfactor = options.get(key + ".bsemin").as<double>();
    bse_vmin =
        num_of_occlevels - Index(bseminfactor * double(num_of_occlevels)) - 1;

    double bsemaxfactor = options.get(key + ".bsemax").as<double>();
    bse_cmax =
        num_of_occlevels + Index(bsemaxfactor * double(num_of_occlevels)) - 1;

  } else if (ranges == "explicit") {
    // get explicit numbers
    rpamax = options.get(key + ".rpamax").as<Index>();
    qpmin = options.get(key + ".qpmin").as<Index>();
    qpmax = options.get(key + ".qpmax").as<Index>();
    bse_vmin = options.get(key + ".bsemin").as<Index>();
    bse_cmax = options.get(key + ".bsemax").as<Index>();
  } else if (ranges == "default") {
    rpamax = num_of_levels - 1;
    qpmin = 0;
    qpmax = 3 * homo + 1;
    bse_vmin = 0;
    bse_cmax = 3 * homo + 1;
  } else if (ranges == "full") {
    rpamax = num_of_levels - 1;
    qpmin = 0;
    qpmax = num_of_levels - 1;
    bse_vmin = 0;
    bse_cmax = num_of_levels - 1;
  }
  std::string ignore_corelevels =
      options.get(key + ".ignore_corelevels").as<std::string>();

  if (ignore_corelevels == "RPA" || ignore_corelevels == "GW" ||
      ignore_corelevels == "BSE") {
    Index ignored_corelevels = CountCoreLevels();
    if (ignore_corelevels == "RPA") {
      rpamin = ignored_corelevels;
    }
    if (ignore_corelevels == "GW" || ignore_corelevels == "RPA") {
      if (qpmin < ignored_corelevels) {
        qpmin = ignored_corelevels;
      }
    }
    if (ignore_corelevels == "GW" || ignore_corelevels == "RPA" ||
        ignore_corelevels == "BSE") {
      if (bse_vmin < ignored_corelevels) {
        bse_vmin = ignored_corelevels;
      }
    }

    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Ignoring " << ignored_corelevels
        << " core levels for " << ignore_corelevels << " and beyond." << flush;
  }

  // check maximum and minimum sizes
  if (rpamax > num_of_levels) {
    rpamax = num_of_levels - 1;
  }
  if (qpmax > num_of_levels) {
    qpmax = num_of_levels - 1;
  }
  if (bse_cmax > num_of_levels) {
    bse_cmax = num_of_levels - 1;
  }
  if (bse_vmin < 0) {
    bse_vmin = 0;
  }
  if (qpmin < 0) {
    qpmin = 0;
  }

  _gwopt.homo = homo;
  _gwopt.qpmin = qpmin;
  _gwopt.qpmax = qpmax;
  _gwopt.rpamin = rpamin;
  _gwopt.rpamax = rpamax;

  _bseopt.vmin = bse_vmin;
  _bseopt.cmax = bse_cmax;
  _bseopt.homo = homo;
  _bseopt.qpmin = qpmin;
  _bseopt.qpmax = qpmax;
  _bseopt.rpamin = rpamin;
  _bseopt.rpamax = rpamax;

  _orbitals.setRPAindices(rpamin, rpamax);
  _orbitals.setGWindices(qpmin, qpmax);
  _orbitals.setBSEindices(bse_vmin, bse_cmax);
  _orbitals.SetFlagUseHqpOffdiag(_bseopt.use_Hqp_offdiag);

  Index bse_vmax = homo;
  Index bse_cmin = homo + 1;
  Index bse_vtotal = bse_vmax - bse_vmin + 1;
  Index bse_ctotal = bse_cmax - bse_cmin + 1;
  Index bse_size = bse_vtotal * bse_ctotal;

  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " RPA level range [" << rpamin
                              << ":" << rpamax << "]" << flush;
  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " GW  level range [" << qpmin
                              << ":" << qpmax << "]" << flush;
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " BSE level range occ[" << bse_vmin << ":" << bse_vmax
      << "]  virt[" << bse_cmin << ":" << bse_cmax << "]" << flush;
  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " BSE Hamiltonian has size "
                              << bse_size << "x" << bse_size << flush;

  _gwopt.reset_3c = options.get(key + ".rebuild_threecenter_freq").as<Index>();

  _bseopt.nmax = options.get(key + ".exctotal").as<Index>();
  if (_bseopt.nmax > bse_size || _bseopt.nmax < 0) {
    _bseopt.nmax = bse_size;
  }

  // eigensolver options
  _bseopt.davidson = options.get(key + ".eigensolver.dodavidson").as<bool>();

  if (_bseopt.davidson) {

    _bseopt.matrixfree =
        options.get(key + ".eigensolver.domatrixfree").as<bool>();

    _bseopt.davidson_correction =
        options.get(key + ".eigensolver.davidson_correction").as<std::string>();

    _bseopt.davidson_ortho =
        options.get(key + ".eigensolver.davidson_ortho").as<std::string>();

    _bseopt.davidson_tolerance =
        options.get(key + ".eigensolver.davidson_tolerance").as<std::string>();

    _bseopt.davidson_update =
        options.get(key + ".eigensolver.davidson_update").as<std::string>();

    _bseopt.davidson_maxiter =
        options.get(key + ".eigensolver.davidson_maxiter").as<Index>();

    // check size
    if (_bseopt.nmax > bse_size / 4) {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp()
          << " Warning : Too many eigenvalues required for Davidson. Default "
             "to Lapack diagonalization"
          << flush;
      _bseopt.davidson = false;
    }
  }

  _bseopt.useTDA = options.get(key + ".useTDA").as<bool>();
  _orbitals.setTDAApprox(_bseopt.useTDA);
  if (!_bseopt.useTDA) {
    XTP_LOG(Log::error, *_pLog) << " BSE type: full" << flush;
  } else {
    XTP_LOG(Log::error, *_pLog) << " BSE type: TDA" << flush;
  }

  _bseopt.use_Hqp_offdiag = options.get(key + ".use_Hqp_offdiag").as<bool>();

  if (!_bseopt.use_Hqp_offdiag) {
    XTP_LOG(Log::error, *_pLog)
        << " BSE without Hqp offdiagonal elements" << flush;
  } else {
    XTP_LOG(Log::error, *_pLog)
        << " BSE with Hqp offdiagonal elements" << flush;
  }

  _bseopt.max_dyn_iter =
      options.get(key + ".dynamical_screening_max_iterations").as<Index>();
  _bseopt.dyn_tolerance =
      options.get(key + ".dynamical_screening_tolerance").as<double>();
  if (_bseopt.max_dyn_iter > 0) {
    _do_dynamical_screening_bse = true;
  }

  _functional = options.get(key + ".vxc.functional").as<std::string>();
  _grid = options.get(key + ".vxc.grid").as<std::string>();

  _auxbasis_name = options.get(key + ".auxbasisset").as<std::string>();
  _dftbasis_name = options.get(key + ".basisset").as<std::string>();
  if (_dftbasis_name != _orbitals.getDFTbasisName()) {
    throw std::runtime_error(
        "Name of the Basisset from .orb file: " + _orbitals.getDFTbasisName() +
        " and from GWBSE optionfile " + _dftbasis_name + " do not agree.");
  }

  std::string mode = options.get(key + ".mode").as<std::string>();
  if (mode == "G0W0") {
    _gwopt.gw_sc_max_iterations = 1;
  } else if (mode == "evGW") {
    _gwopt.g_sc_limit = 0.1 * _gwopt.gw_sc_limit;
    _gwopt.eta = 0.1;
  }

  // Check if the exact exchange of the functional matches the exact exchange in
  // orb file. If not display a warning.
  double ScaHFX_temp = Vxc_Potential<Vxc_Grid>::getExactExchange(_functional);
  if (ScaHFX_temp != _orbitals.getScaHFX()) {
    XTP_LOG(Log::error, *_pLog)
        << (boost::format(
                "WARNING: GWBSE exact exchange a=%s differs from qmpackage "
                "exact exchange a=%s, \n probably your functionals are "
                "inconsistent \n or dft orbitals were loaded from a package"
                " other than votca or orca. \n The GWBSE exact exchange will "
                "be used.") %
            ScaHFX_temp % _orbitals.getScaHFX())
               .str()
        << std::flush;
    _orbitals.setScaHFX(ScaHFX_temp);
  }

  XTP_LOG(Log::error, *_pLog) << " Running GW as: " << mode << flush;
  _gwopt.ScaHFX = _orbitals.getScaHFX();

  _gwopt.shift = options.get(key + ".scissor_shift").as<double>();
  _gwopt.g_sc_limit =
      options.get(key + ".g_sc_limit").as<double>();  // convergence criteria
                                                      // for qp iteration
                                                      // [Hartree]]
  _gwopt.g_sc_max_iterations =
      options.get(key + ".g_sc_max_iterations").as<Index>();  // convergence
                                                              // criteria for qp
                                                              // iteration
                                                              // [Hartree]]

  if (mode == "evGW") {
    _gwopt.gw_sc_max_iterations =
        options.get(key + ".gw_sc_max_iterations").as<Index>();
  }

  _gwopt.gw_sc_limit =
      options.get(key + ".gw_sc_limit").as<double>();  // convergence criteria
                                                       // for shift it
  XTP_LOG(Log::error, *_pLog)
      << " g_sc_limit [Hartree]: " << _gwopt.g_sc_limit << flush;
  if (_gwopt.gw_sc_max_iterations > 1) {
    XTP_LOG(Log::error, *_pLog)
        << " gw_sc_limit [Hartree]: " << _gwopt.gw_sc_limit << flush;
  }
  _bseopt.min_print_weight =
      options.get(key + ".bse_print_weight").as<double>();
  // print exciton WF composition weight larger than minimum

  // possible tasks
  std::string tasks_string = options.get(key + ".tasks").as<std::string>();
  boost::algorithm::to_lower(tasks_string);
  if (tasks_string.find("all") != std::string::npos) {
    _do_gw = true;
    _do_bse_singlets = true;
    _do_bse_triplets = true;
  }
  if (tasks_string.find("gw") != std::string::npos) {
    _do_gw = true;
  }
  if (tasks_string.find("singlets") != std::string::npos) {
    _do_bse_singlets = true;
  }
  if (tasks_string.find("triplets") != std::string::npos) {
    _do_bse_triplets = true;
  }

  XTP_LOG(Log::error, *_pLog) << " Tasks: " << flush;
  if (_do_gw) {
    XTP_LOG(Log::error, *_pLog) << " GW " << flush;
  }
  if (_do_bse_singlets) {
    XTP_LOG(Log::error, *_pLog) << " singlets " << flush;
  }
  if (_do_bse_triplets) {
    XTP_LOG(Log::error, *_pLog) << " triplets " << flush;
  }
  XTP_LOG(Log::error, *_pLog) << " Store: " << flush;
  if (_do_gw) {
    XTP_LOG(Log::error, *_pLog) << " GW " << flush;
  }

  if (options.exists(key + ".fragments")) {
    std::vector<tools::Property*> prop_region =
        options.Select(key + ".fragments.fragment");
    Index index = 0;
    for (tools::Property* prop : prop_region) {
      std::string indices =
          prop->ifExistsReturnElseThrowRuntimeError<std::string>("indices");
      _fragments.push_back(QMFragment<BSE_Population>(index, indices));
      index++;
    }
  }

  _gwopt.sigma_integration =
      options.get(key + ".sigma_integrator").as<std::string>();
  XTP_LOG(Log::error, *_pLog)
      << " Sigma integration: " << _gwopt.sigma_integration << flush;
  _gwopt.eta = options.get(key + ".eta").as<double>();
  XTP_LOG(Log::error, *_pLog) << " eta: " << _gwopt.eta << flush;
  if (_gwopt.sigma_integration == "exact") {
    XTP_LOG(Log::error, *_pLog)
        << " RPA Hamiltonian size: " << (homo + 1 - rpamin) * (rpamax - homo)
        << flush;
  }
  _gwopt.order = options.get(key + ".quadrature_order").as<Index>();
  XTP_LOG(Log::error, *_pLog)
      << " Quadrature integration order : " << _gwopt.order << flush;
  _gwopt.quadrature_scheme =
      options.get(key + ".quadrature_scheme").as<std::string>();
  XTP_LOG(Log::error, *_pLog)
      << " Quadrature integration scheme : " << _gwopt.quadrature_scheme
      << flush;
  _gwopt.alpha = options.get(key + ".alpha").as<double>();
  XTP_LOG(Log::error, *_pLog)
      << " Alpha smoothing parameter : " << _gwopt.alpha << flush;

  _gwopt.qp_solver = options.get(key + ".qp_solver").as<std::string>();
  _gwopt.qp_grid_steps = options.get(key + ".qp_grid_steps").as<Index>();
  _gwopt.qp_grid_spacing = options.get(key + ".qp_grid_spacing").as<double>();
  XTP_LOG(Log::error, *_pLog) << " QP solver: " << _gwopt.qp_solver << flush;
  if (_gwopt.qp_solver == "grid") {
    XTP_LOG(Log::error, *_pLog)
        << " QP grid steps: " << _gwopt.qp_grid_steps << flush;
    XTP_LOG(Log::error, *_pLog)
        << " QP grid spacing: " << _gwopt.qp_grid_spacing << flush;
  }
  _gwopt.gw_mixing_order =
      options.get(key + ".gw_mixing_order").as<Index>();  // max history in
                                                          // mixing (0: plain,
                                                          // 1: linear, >1
                                                          // Anderson)

  _gwopt.gw_mixing_alpha = options.get(key + ".gw_mixing_alpha").as<double>();

  if (mode == "evGW") {
    if (_gwopt.gw_mixing_order == 0) {
      XTP_LOG(Log::error, *_pLog) << " evGW with plain update " << std::flush;
    } else if (_gwopt.gw_mixing_order == 1) {
      XTP_LOG(Log::error, *_pLog) << " evGW with linear update using alpha "
                                  << _gwopt.gw_mixing_alpha << std::flush;
    } else {
      XTP_LOG(Log::error, *_pLog) << " evGW with Anderson update with history "
                                  << _gwopt.gw_mixing_order << " using alpha "
                                  << _gwopt.gw_mixing_alpha << std::flush;
    }
  }

  _sigma_plot_states =
      options.get(key + ".sigma_plot_states").as<std::string>();
  _sigma_plot_steps = options.get(key + ".sigma_plot_steps").as<Index>();
  _sigma_plot_spacing = options.get(key + ".sigma_plot_spacing").as<double>();
  _sigma_plot_filename =
      options.get(key + ".sigma_plot_filename").as<std::string>();
  if (!_sigma_plot_states.empty()) {
    XTP_LOG(Log::error, *_pLog)
        << " Sigma plot states: " << _sigma_plot_states << flush;
    XTP_LOG(Log::error, *_pLog)
        << " Sigma plot steps: " << _sigma_plot_steps << flush;
    XTP_LOG(Log::error, *_pLog)
        << " Sigma plot spacing: " << _sigma_plot_spacing << flush;
    XTP_LOG(Log::error, *_pLog)
        << " Sigma plot filename: " << _sigma_plot_filename << flush;
  }
}

void GWBSE::addoutput(tools::Property& summary) {

  const double hrt2ev = tools::conv::hrt2ev;
  tools::Property& gwbse_summary = summary.add("GWBSE", "");
  if (_do_gw) {
    gwbse_summary.setAttribute("units", "eV");
    gwbse_summary.setAttribute(
        "DFTEnergy",
        (format("%1$+1.6f ") % (_orbitals.getDFTTotalEnergy() * hrt2ev)).str());

    tools::Property& dft_summary = gwbse_summary.add("dft", "");
    dft_summary.setAttribute("HOMO", _gwopt.homo);
    dft_summary.setAttribute("LUMO", _gwopt.homo + 1);

    for (Index state = 0; state < _gwopt.qpmax + 1 - _gwopt.qpmin; state++) {
      tools::Property& level_summary = dft_summary.add("level", "");
      level_summary.setAttribute("number", state + _gwopt.qpmin);
      level_summary.add(
          "dft_energy",
          (format("%1$+1.6f ") %
           (_orbitals.MOs().eigenvalues()(state + _gwopt.qpmin) * hrt2ev))
              .str());
      level_summary.add(
          "gw_energy",
          (format("%1$+1.6f ") % (_orbitals.QPpertEnergies()(state) * hrt2ev))
              .str());

      level_summary.add("qp_energy",
                        (format("%1$+1.6f ") %
                         (_orbitals.QPdiag().eigenvalues()(state) * hrt2ev))
                            .str());
    }
  }
  if (_do_bse_singlets) {
    tools::Property& singlet_summary = gwbse_summary.add("singlets", "");
    for (Index state = 0; state < _bseopt.nmax; ++state) {
      tools::Property& level_summary = singlet_summary.add("level", "");
      level_summary.setAttribute("number", state + 1);
      level_summary.add(
          "omega", (format("%1$+1.6f ") %
                    (_orbitals.BSESinglets().eigenvalues()(state) * hrt2ev))
                       .str());
      if (_orbitals.hasTransitionDipoles()) {

        const Eigen::Vector3d& dipoles = (_orbitals.TransitionDipoles())[state];
        double f = 2 * dipoles.squaredNorm() *
                   _orbitals.BSESinglets().eigenvalues()(state) / 3.0;

        level_summary.add("f", (format("%1$+1.6f ") % f).str());
        tools::Property& dipol_summary = level_summary.add(
            "Trdipole", (format("%1$+1.4f %2$+1.4f %3$+1.4f") % dipoles.x() %
                         dipoles.y() % dipoles.z())
                            .str());
        dipol_summary.setAttribute("unit", "e*bohr");
        dipol_summary.setAttribute("gauge", "length");
      }
    }
  }
  if (_do_bse_triplets) {
    tools::Property& triplet_summary = gwbse_summary.add("triplets", "");
    for (Index state = 0; state < _bseopt.nmax; ++state) {
      tools::Property& level_summary = triplet_summary.add("level", "");
      level_summary.setAttribute("number", state + 1);
      level_summary.add(
          "omega", (format("%1$+1.6f ") %
                    (_orbitals.BSETriplets().eigenvalues()(state) * hrt2ev))
                       .str());
    }
  }
  return;
}

/*
 *    Many-body Green's fuctions theory implementation
 *
 *  data required from orbitals file
 *  - atomic coordinates
 *  - DFT molecular orbitals (energies and coeffcients)
 *  - DFT exchange-correlation potential matrix in atomic orbitals
 *  - number of electrons, number of levels
 */

Eigen::MatrixXd GWBSE::CalculateVXC(const AOBasis& dftbasis) {
  if (_orbitals.getXCFunctionalName().empty()) {
    _orbitals.setXCFunctionalName(_functional);
  } else {
    if (!(_functional == _orbitals.getXCFunctionalName())) {
      throw std::runtime_error("Functionals from DFT " +
                               _orbitals.getXCFunctionalName() + " GWBSE " +
                               _functional + " differ!");
    }
  }

  Vxc_Grid grid;
  grid.GridSetup(_grid, _orbitals.QMAtoms(), dftbasis);
  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Setup grid for integration with gridsize: " << _grid
      << " with " << grid.getGridSize() << " points, divided into "
      << grid.getBoxesSize() << " boxes" << flush;
  Vxc_Potential<Vxc_Grid> vxcpotential(grid);
  vxcpotential.setXCfunctional(_functional);
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Integrating Vxc in VOTCA with functional "
      << _functional << flush;
  Eigen::MatrixXd DMAT = _orbitals.DensityMatrixGroundState();
  Mat_p_Energy e_vxc_ao = vxcpotential.IntegrateVXC(DMAT);
  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Calculated Vxc in VOTCA" << flush;
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Set hybrid exchange factor: " << _orbitals.getScaHFX()
      << flush;
  Index qptotal = _gwopt.qpmax - _gwopt.qpmin + 1;
  Index basissize = Index(_orbitals.MOs().eigenvectors().rows());
  Eigen::MatrixXd mos =
      _orbitals.MOs().eigenvectors().block(0, _gwopt.qpmin, basissize, qptotal);

  Eigen::MatrixXd vxc = mos.transpose() * e_vxc_ao.matrix() * mos;
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Calculated exchange-correlation expectation values "
      << flush;

  return vxc;
}

bool GWBSE::Evaluate() {

  // set the parallelization
  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Using "
                              << OPENMP::getMaxThreads() << " threads" << flush;

  if (XTP_HAS_MKL_OVERLOAD()) {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Using MKL overload for Eigen " << flush;
  } else {
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp()
        << " Using native Eigen implementation, no BLAS overload " << flush;
  }

  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Molecule Coordinates [A] " << flush;
  for (QMAtom& atom : _orbitals.QMAtoms()) {
    std::string output =
        (boost::format("  %1$s"
                       "   %2$+1.4f %3$+1.4f %4$+1.4f") %
         atom.getElement() % (atom.getPos().x() * tools::conv::bohr2ang) %
         (atom.getPos().y() * tools::conv::bohr2ang) %
         (atom.getPos().z() * tools::conv::bohr2ang))
            .str();

    XTP_LOG(Log::error, *_pLog) << output << flush;
  }

  std::string dft_package = _orbitals.getQMpackage();
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " DFT data was created by " << dft_package << flush;

  BasisSet dftbs;
  dftbs.Load(_dftbasis_name);

  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Loaded DFT Basis Set " << _dftbasis_name << flush;

  // fill DFT AO basis by going through all atoms
  AOBasis dftbasis;
  dftbasis.Fill(dftbs, _orbitals.QMAtoms());
  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Filled DFT Basis of size "
                              << dftbasis.AOBasisSize() << flush;

  // load auxiliary basis set (element-wise information) from xml file
  BasisSet auxbs;
  auxbs.Load(_auxbasis_name);
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Loaded Auxbasis Set " << _auxbasis_name << flush;

  // fill auxiliary AO basis by going through all atoms
  AOBasis auxbasis;
  auxbasis.Fill(auxbs, _orbitals.QMAtoms());
  _orbitals.setAuxbasisName(_auxbasis_name);
  XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Filled Auxbasis of size "
                              << auxbasis.AOBasisSize() << flush;

  if ((_do_bse_singlets || _do_bse_triplets) && _fragments.size() > 0) {
    for (const auto& frag : _fragments) {
      XTP_LOG(Log::error, *_pLog) << TimeStamp() << " Fragment " << frag.getId()
                                  << " size:" << frag.size() << flush;
    }
  }

  if (!_do_gw && !_orbitals.hasQPdiag()) {
    throw std::runtime_error(
        "You want no GW calculation but the orb file has no QPcoefficients for "
        "BSE");
  }
  TCMatrix_gwbse Mmn(*_pLog);
  // rpamin here, because RPA needs till rpamin
  Index max_3c = std::max(_bseopt.cmax, _gwopt.qpmax);
  Mmn.Initialize(auxbasis.AOBasisSize(), _gwopt.rpamin, max_3c, _gwopt.rpamin,
                 _gwopt.rpamax);
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp()
      << " Calculating Mmn_beta (3-center-repulsion x orbitals)  " << flush;
  Mmn.Fill(auxbasis, dftbasis, _orbitals.MOs().eigenvectors());
  XTP_LOG(Log::info, *_pLog)
      << TimeStamp() << " Removed " << Mmn.Removedfunctions()
      << " functions from Aux Coulomb matrix to avoid near linear dependencies"
      << flush;
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " Calculated Mmn_beta (3-center-repulsion x orbitals)  "
      << flush;

  Eigen::MatrixXd Hqp;
  if (_do_gw) {
    Eigen::MatrixXd vxc = CalculateVXC(dftbasis);
    GW gw = GW(*_pLog, Mmn, vxc, _orbitals.MOs().eigenvalues());
    gw.configure(_gwopt);
    gw.CalculateGWPerturbation();

    if (!_sigma_plot_states.empty()) {
      gw.PlotSigma(_sigma_plot_filename, _sigma_plot_steps, _sigma_plot_spacing,
                   _sigma_plot_states);
    }

    // store perturbative QP energy data in orbitals object (DFT, S_x,S_c, V_xc,
    // E_qp)
    _orbitals.QPpertEnergies() = gw.getGWAResults();
    _orbitals.RPAInputEnergies() = gw.RPAInputEnergies();

    XTP_LOG(Log::info, *_pLog)
        << TimeStamp() << " Calculating offdiagonal part of Sigma  " << flush;
    gw.CalculateHQP();
    XTP_LOG(Log::error, *_pLog)
        << TimeStamp() << " Calculated offdiagonal part of Sigma  " << flush;

    Hqp = gw.getHQP();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es =
        gw.DiagonalizeQPHamiltonian();
    if (es.info() == Eigen::ComputationInfo::Success) {
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Diagonalized QP Hamiltonian  " << flush;
    }

    _orbitals.QPdiag().eigenvectors() = es.eigenvectors();
    _orbitals.QPdiag().eigenvalues() = es.eigenvalues();
  } else {
    if (_orbitals.getGWAmax() != _gwopt.qpmax ||
        _orbitals.getGWAmin() != _gwopt.qpmin ||
        _orbitals.getRPAmax() != _gwopt.rpamax ||
        _orbitals.getRPAmin() != _gwopt.rpamin) {
      throw std::runtime_error(
          "The ranges for GW and RPA do not agree with the ranges from the "
          ".orb file, rerun your GW calculation");
    }
    const Eigen::MatrixXd& qpcoeff = _orbitals.QPdiag().eigenvectors();

    Hqp = qpcoeff * _orbitals.QPdiag().eigenvalues().asDiagonal() *
          qpcoeff.transpose();
  }

  // proceed only if BSE requested
  if (_do_bse_singlets || _do_bse_triplets) {

    BSE bse = BSE(*_pLog, Mmn);
    bse.configure(_bseopt, _orbitals.RPAInputEnergies(), Hqp);

    // store the direct contribution to the static BSE results
    Eigen::VectorXd Hd_static_contrib_triplet;
    Eigen::VectorXd Hd_static_contrib_singlet;

    if (_do_bse_triplets) {
      bse.Solve_triplets(_orbitals);
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Solved BSE for triplets " << flush;
      Hd_static_contrib_triplet = bse.Analyze_triplets(_fragments, _orbitals);
    }

    if (_do_bse_singlets) {
      bse.Solve_singlets(_orbitals);
      XTP_LOG(Log::error, *_pLog)
          << TimeStamp() << " Solved BSE for singlets " << flush;
      Hd_static_contrib_singlet = bse.Analyze_singlets(_fragments, _orbitals);
    }

    // do perturbative dynamical screening in BSE
    if (_do_dynamical_screening_bse) {

      if (_do_bse_triplets) {
        bse.Perturbative_DynamicalScreening(QMStateType(QMStateType::Triplet),
                                            _orbitals,
                                            Hd_static_contrib_triplet);
      }

      if (_do_bse_singlets) {
        bse.Perturbative_DynamicalScreening(QMStateType(QMStateType::Singlet),
                                            _orbitals,
                                            Hd_static_contrib_singlet);
      }
    }
  }
  XTP_LOG(Log::error, *_pLog)
      << TimeStamp() << " GWBSE calculation finished " << flush;
  return true;
}

}  // namespace xtp
}  // namespace votca
