/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE gw_test

// Standard includes
#include <algorithm>
#include <fstream>
#include <sstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/environmentscreening.h"
#include "votca/xtp/gw.h"
#include "votca/xtp/gw_uks.h"
#include "votca/xtp/ppm.h"
#include "votca/xtp/rpa_uks.h"
#include "votca/xtp/sigmafactory_uks.h"

// for SetSharedPPM
#include "../libxtp/self_energy_evaluators/sigma_ppm_uks.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigma_base.h"
#include "votca/xtp/sigmafactory.h"

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;
using namespace std;

namespace {

struct GWTestSystem {
  Eigen::VectorXd mo_eigenvalues;
  Eigen::MatrixXd mo_eigenvectors;
  Eigen::MatrixXd vxc;
  Logger log;
  TCMatrix_gwbse Mmn;

  GWTestSystem(const std::string& mo_file, const std::string& vxc_file)
      : mo_eigenvalues(Eigen::VectorXd::Zero(17)) {
    mo_eigenvalues << -10.6784, -0.746424, -0.394948, -0.394948, -0.394948,
        0.165212, 0.227713, 0.227713, 0.227713, 0.763971, 0.763971, 0.763971,
        1.05054, 1.13372, 1.13372, 1.13372, 1.72964;

    mo_eigenvectors = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/gw/" + mo_file);
    vxc = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/gw/" + vxc_file);

    Orbitals orbitals;
    orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                    "/gw/molecule.xyz");

    BasisSet basis;
    basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
    orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");

    AOBasis aobasis;
    aobasis.Fill(basis, orbitals.QMAtoms());

    orbitals.setNumberOfOccupiedLevels(4);
    orbitals.MOs().eigenvectors() = mo_eigenvectors;

    Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
    Mmn.Fill(aobasis, aobasis, mo_eigenvectors);
  }
};

GW::options MakeGWTestOptions() {
  GW::options opt;
  opt.ScaHFX = 0;
  opt.homo = 4;
  opt.qpmax = 16;
  opt.qpmin = 0;
  opt.rpamax = 16;
  opt.rpamin = 0;
  opt.gw_sc_max_iterations = 1;
  opt.qp_solver = "grid";
  opt.eta = 1e-3;
  opt.sigma_integration = "ppm";
  opt.reset_3c = 5;
  opt.gw_mixing_order = 0;
  opt.gw_mixing_alpha = 0.7;
  opt.g_sc_limit = 1e-5;
  opt.g_sc_max_iterations = 50;
  opt.gw_sc_limit = 1e-5;
  // the tests below were written for these; the defaults are adaptive and
  // continuity (tests of those set them explicitly)
  opt.screening_update = "every";
  opt.qp_root_continuity = false;

  return opt;
}

Eigen::VectorXd RunGWPerturbation(const GWTestSystem& system,
                                  const GW::options& opt) {
  GW gw(const_cast<Logger&>(system.log),
        const_cast<TCMatrix_gwbse&>(system.Mmn), system.vxc,
        system.mo_eigenvalues);
  gw.configure(opt);
  gw.CalculateGWPerturbation();
  return gw.getGWAResults();
}

}  // namespace

BOOST_AUTO_TEST_SUITE(gw_test)

BOOST_AUTO_TEST_CASE(gw_full) {
  libint2::initialize();
  Eigen::VectorXd mo_eigenvalues = Eigen::VectorXd::Zero(17);
  mo_eigenvalues << -10.6784, -0.746424, -0.394948, -0.394948, -0.394948,
      0.165212, 0.227713, 0.227713, 0.227713, 0.763971, 0.763971, 0.763971,
      1.05054, 1.13372, 1.13372, 1.13372, 1.72964;
  Eigen::MatrixXd mo_eigenvectors =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/gw/mo_eigenvectors.mm");

  Eigen::MatrixXd vxc = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gw/vxc.mm");

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/gw/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  orbitals.setNumberOfOccupiedLevels(4);
  orbitals.MOs().eigenvectors() = mo_eigenvectors;
  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, mo_eigenvectors);

  GW::options opt;
  opt.ScaHFX = 0;
  opt.homo = 4;
  opt.qpmax = 16;
  opt.qpmin = 0;
  opt.rpamax = 16;
  opt.rpamin = 0;
  opt.gw_sc_max_iterations = 1;
  opt.eta = 1e-3;
  opt.sigma_integration = "ppm";
  opt.reset_3c = 5;
  opt.qp_solver = "grid";
  opt.qp_grid_steps = 601;
  opt.qp_grid_spacing = 0.005;
  opt.gw_mixing_order = 0;
  opt.gw_mixing_alpha = 0.7;
  opt.g_sc_limit = 1e-5;
  opt.g_sc_max_iterations = 50;
  opt.gw_sc_limit = 1e-5;

  GW gw(log, Mmn, vxc, mo_eigenvalues);
  gw.configure(opt);

  Eigen::MatrixXd ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gw/ref.mm");

  gw.CalculateGWPerturbation();
  Eigen::VectorXd diag = gw.getGWAResults();
  bool check_diag = ref.diagonal().isApprox(diag, 1e-4);
  if (!check_diag) {
    cout << "GW energies" << endl;
    cout << diag << endl;
    cout << "GW energies ref" << endl;
    cout << ref.diagonal() << endl;
  }
  BOOST_CHECK_EQUAL(check_diag, true);

  gw.CalculateHQP();
  Eigen::MatrixXd offdiag = gw.getHQP();

  bool check_offdiag = ref.isApprox(offdiag, 1e-4);
  if (!check_offdiag) {
    cout << "GW energies" << endl;
    cout << offdiag << endl;
    cout << "GW energies ref" << endl;
    cout << ref << endl;
  }
  BOOST_CHECK_EQUAL(check_offdiag, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(gw_full_QP_grid) {
  libint2::initialize();
  Eigen::VectorXd mo_eigenvalues = Eigen::VectorXd::Zero(17);
  mo_eigenvalues << -10.6784, -0.746424, -0.394948, -0.394948, -0.394948,
      0.165212, 0.227713, 0.227713, 0.227713, 0.763971, 0.763971, 0.763971,
      1.05054, 1.13372, 1.13372, 1.13372, 1.72964;
  Eigen::MatrixXd mo_eigenvectors =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/gw/mo_eigenvectors2.mm");

  Eigen::MatrixXd vxc = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gw/vxc2.mm");

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/gw/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  orbitals.setNumberOfOccupiedLevels(4);
  orbitals.MOs().eigenvectors() = mo_eigenvectors;
  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, mo_eigenvectors);

  GW::options opt;
  opt.ScaHFX = 0;
  opt.homo = 4;
  opt.qpmax = 16;
  opt.qpmin = 0;
  opt.rpamax = 16;
  opt.rpamin = 0;
  opt.gw_sc_max_iterations = 1;
  opt.qp_solver = "grid";
  opt.eta = 1e-3;
  opt.sigma_integration = "ppm";
  opt.reset_3c = 5;
  opt.qp_grid_steps = 601;
  opt.qp_grid_spacing = 0.005;
  opt.gw_mixing_order = 0;
  opt.gw_mixing_alpha = 0.7;
  opt.g_sc_limit = 1e-5;
  opt.g_sc_max_iterations = 50;
  opt.gw_sc_limit = 1e-5;

  GW gw(log, Mmn, vxc, mo_eigenvalues);
  gw.configure(opt);

  gw.CalculateGWPerturbation();
  Eigen::VectorXd diag = gw.getGWAResults();

  Eigen::MatrixXd ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gw/ref2.mm");

  bool check_diag = ref.diagonal().isApprox(diag, 1e-4);
  if (!check_diag) {
    cout << "GW energies" << endl;
    cout << diag << endl;
    cout << "GW energies ref" << endl;
    cout << ref.diagonal() << endl;
  }
  BOOST_CHECK_EQUAL(check_diag, true);

  gw.CalculateHQP();
  Eigen::MatrixXd offdiag = gw.getHQP();

  bool check_offdiag = ref.isApprox(offdiag, 1e-4);
  if (!check_offdiag) {
    cout << "GW energies" << endl;
    cout << offdiag << endl;
    cout << "GW energies ref" << endl;
    cout << ref << endl;
  }
  BOOST_CHECK_EQUAL(check_offdiag, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(gw_full_QP_grid_canonical_options) {
  libint2::initialize();

  GWTestSystem system("mo_eigenvectors2.mm", "vxc2.mm");

  GW::options opt = MakeGWTestOptions();

  // New canonical QP-search controls. These reproduce the legacy dense window:
  // old qp_grid_steps = 601 and qp_grid_spacing = 0.005 imply
  // full half-width = 0.5 * (601 - 1) * 0.005 = 1.5 Ha
  // dense spacing    = 0.005 Ha
  // legacy adaptive shell width would be 3.0 / (150 - 1), but for this test
  // we intentionally use the canonical controls directly.
  opt.qp_full_window_half_width = 1.5;
  opt.qp_dense_spacing = 0.005;
  opt.qp_adaptive_shell_width = 0.02;
  opt.qp_adaptive_shell_count = 0;
  opt.qp_grid_search_mode = "adaptive_with_dense_fallback";
  opt.qp_root_finder = "bisection";

  Eigen::VectorXd diag = RunGWPerturbation(system, opt);

  Eigen::MatrixXd ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gw/ref2.mm");

  bool check_diag = ref.diagonal().isApprox(diag, 1e-4);
  if (!check_diag) {
    cout << "GW energies (canonical controls)" << endl;
    cout << diag << endl;
    cout << "GW energies ref" << endl;
    cout << ref.diagonal() << endl;
  }
  BOOST_CHECK_EQUAL(check_diag, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(gw_full_QP_grid_brent_matches_bisection) {
  libint2::initialize();

  GWTestSystem system("mo_eigenvectors2.mm", "vxc2.mm");

  GW::options opt_bisect = MakeGWTestOptions();
  opt_bisect.qp_full_window_half_width = 1.5;
  opt_bisect.qp_dense_spacing = 0.005;
  opt_bisect.qp_adaptive_shell_width = 0.02;
  opt_bisect.qp_adaptive_shell_count = 0;
  opt_bisect.qp_grid_search_mode = "adaptive_with_dense_fallback";
  opt_bisect.qp_root_finder = "bisection";

  GW::options opt_brent = opt_bisect;
  opt_brent.qp_root_finder = "brent";

  Eigen::VectorXd diag_bisect = RunGWPerturbation(system, opt_bisect);
  Eigen::VectorXd diag_brent = RunGWPerturbation(system, opt_brent);

  bool check_diag = diag_bisect.isApprox(diag_brent, 1e-5);
  if (!check_diag) {
    cout << "GW energies (bisection)" << endl;
    cout << diag_bisect << endl;
    cout << "GW energies (brent)" << endl;
    cout << diag_brent << endl;
    cout << "Difference" << endl;
    cout << (diag_bisect - diag_brent) << endl;
  }
  BOOST_CHECK_EQUAL(check_diag, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(qsgw_ppm) {
  // QSGW self-consistency with PPM sigma on methane/3-21G.
  // Uses qpmax=16 with qsgw_max_virt_correction=0.35 Ha to auto-trim level 16.
  // Checks: convergence, starting-point independence, rotation unitarity.
  // Note: libint2::initialize() is not called here because some platforms
  // crash on repeated init/finalize cycles. libint2 is already initialized
  // by the preceding test suite (gw_full initializes it at suite start).
  // We call initialize only if not already done.
  if (!libint2::initialized()) libint2::initialize();

  Eigen::VectorXd mo_eigenvalues = Eigen::VectorXd::Zero(17);
  mo_eigenvalues << -10.6784, -0.746424, -0.394948, -0.394948, -0.394948,
      0.165212, 0.227713, 0.227713, 0.227713, 0.763971, 0.763971, 0.763971,
      1.05054, 1.13372, 1.13372, 1.13372, 1.72964;
  Eigen::MatrixXd mo_eigenvectors =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/gw/mo_eigenvectors.mm");
  Eigen::MatrixXd vxc = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gw/vxc.mm");
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/gw/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  GW::options opt = MakeGWTestOptions();
  opt.do_qsgw = true;
  opt.qsgw_max_iterations = 50;
  opt.qsgw_sc_limit = 1e-5;
  opt.gw_mixing_order = 20;
  opt.gw_mixing_alpha = 0.2;
  opt.qpmax = 13;
  opt.qsgw_max_virt_correction = 0.35;

  // G0W0 seed
  opt.gw_sc_max_iterations = 1;
  Logger log1;
  TCMatrix_gwbse Mmn1;
  Mmn1.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn1.Fill(aobasis, aobasis, mo_eigenvectors);
  // Slice vxc to QP window — gwbse.cc does this before constructing GW
  const votca::Index qptotal_gw = opt.qpmax - opt.qpmin + 1;
  Eigen::MatrixXd vxc_gw =
      vxc.block(opt.qpmin, opt.qpmin, qptotal_gw, qptotal_gw);
  GW gw1(log1, Mmn1, vxc_gw, mo_eigenvalues);
  gw1.configure(opt);
  gw1.CalculateGWPerturbation();
  Mmn1.Fill(aobasis, aobasis, mo_eigenvectors);
  gw1.CalculateQSGW();
  Eigen::VectorXd e_g0w0 = gw1.getGWAResults();

  // evGW seed — independent Mmn to avoid shared state
  opt.gw_sc_max_iterations = 50;
  Logger log2;
  TCMatrix_gwbse Mmn2;
  Mmn2.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn2.Fill(aobasis, aobasis, mo_eigenvectors);
  GW gw2(log2, Mmn2, vxc_gw, mo_eigenvalues);
  gw2.configure(opt);
  gw2.CalculateGWPerturbation();
  Mmn2.Fill(aobasis, aobasis, mo_eigenvectors);
  gw2.CalculateQSGW();
  Eigen::VectorXd e_evgw = gw2.getGWAResults();

  // Starting-point independence
  bool check_spi = e_g0w0.isApprox(e_evgw, 1e-4);
  if (!check_spi) {
    cout << "QSGW G0W0: " << e_g0w0.transpose() << endl;
    cout << "QSGW evGW: " << e_evgw.transpose() << endl;
  }
  BOOST_CHECK_EQUAL(check_spi, true);

  // Physical sanity: HOMO negative, LUMO positive, gap positive
  BOOST_CHECK_LT(e_g0w0(4), 0.0);
  BOOST_CHECK_GT(e_g0w0(5), 0.0);
  BOOST_CHECK_GT(e_g0w0(5) - e_g0w0(4), 0.0);

  // Rotation matrix is 17x17 and unitary
  const Eigen::MatrixXd& U = gw1.getQSGWRotation();
  BOOST_CHECK_EQUAL(U.rows(), 14);
  BOOST_CHECK_EQUAL(U.cols(), 14);
  bool check_unitary =
      (U.transpose() * U).isApprox(Eigen::MatrixXd::Identity(14, 14), 1e-6);
  BOOST_CHECK_EQUAL(check_unitary, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(qsgw_virtual_threshold) {
  // Virtual-level threshold auto-trims level 16 and preserves its
  // perturbative seed energy exactly in getGWAResults().
  if (!libint2::initialized()) libint2::initialize();

  Eigen::VectorXd mo_eigenvalues = Eigen::VectorXd::Zero(17);
  mo_eigenvalues << -10.6784, -0.746424, -0.394948, -0.394948, -0.394948,
      0.165212, 0.227713, 0.227713, 0.227713, 0.763971, 0.763971, 0.763971,
      1.05054, 1.13372, 1.13372, 1.13372, 1.72964;
  Eigen::MatrixXd mo_eigenvectors =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/gw/mo_eigenvectors.mm");
  Eigen::MatrixXd vxc = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gw/vxc.mm");
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/gw/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, mo_eigenvectors);

  GW::options opt = MakeGWTestOptions();
  opt.do_qsgw = true;
  opt.qsgw_max_iterations = 50;
  opt.qsgw_sc_limit = 1e-5;
  opt.gw_mixing_order = 20;
  opt.gw_mixing_alpha = 0.2;
  opt.gw_sc_max_iterations = 1;
  opt.qpmax = 16;
  opt.qsgw_max_virt_correction = 0.2;

  // Slice vxc to QP window
  const votca::Index qptotal_vt = opt.qpmax - opt.qpmin + 1;
  Eigen::MatrixXd vxc_vt =
      vxc.block(opt.qpmin, opt.qpmin, qptotal_vt, qptotal_vt);
  GW gw(log, Mmn, vxc_vt, mo_eigenvalues);
  gw.configure(opt);
  gw.CalculateGWPerturbation();
  const double e16_seed = gw.getGWAResults()(16);

  // Rebuild Mmn_ in the DFT-MO basis before QSGW — CalculateGWPerturbation
  // calls PrepareScreening which transforms Mmn_ into the PPM eigenbasis via
  // MultiplyRightWithAuxMatrix. CalculateQSGW must start from the DFT-MO
  // basis (gwbse.cc rebuilds Mmn_ before calling CalculateQSGW in real runs).
  Mmn.Fill(aobasis, aobasis, mo_eigenvectors);
  gw.CalculateQSGW();
  Eigen::VectorXd e_qsgw = gw.getGWAResults();

  // Excluded level 16 retains seed energy exactly
  BOOST_CHECK_CLOSE(e_qsgw(16), e16_seed, 1e-6);

  // Rotation is 17x17 with identity block for excluded level
  const Eigen::MatrixXd& U = gw.getQSGWRotation();
  BOOST_CHECK_EQUAL(U.rows(), 17);
  BOOST_CHECK_EQUAL(U.cols(), 17);
  BOOST_CHECK_CLOSE(U(16, 16), 1.0, 1e-6);

  const Eigen::VectorXd& seed_energies = gw.getQSGWSeedEnergies();
  BOOST_CHECK_EQUAL(seed_energies.size(), 17);

  libint2::finalize();
}

namespace {
// A symmetric, negative semidefinite stand-in for an environment reaction
// field in the auxiliary metric: what EnvironmentScreening produces, minus
// the physics, which the tests below do not need.
Eigen::MatrixXd SomeReactionField(votca::Index n) {
  std::srand(42);
  const Eigen::MatrixXd G = Eigen::MatrixXd::Random(n, n);
  return -0.05 * G * G.transpose() / double(n);
}
}  // namespace

// Step 4. The static COH+SEX self-energy of an environment reaction field,
// against the sum it is defined by, written out term by term:
//   Sigma_nn' = 1/2 sum_m s_m M_nm R M_n'm^T,  s_m = -1 occupied, +1 virtual.
// The reference is taken while Mmn is still in its fill-time frame; GW then
// computes it after nothing has rotated it yet.
BOOST_AUTO_TEST_CASE(reaction_field_self_energy_matches_explicit_sum) {
  if (!libint2::initialized()) libint2::initialize();
  GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
  const GW::options opt = MakeGWTestOptions();
  const Eigen::MatrixXd R = SomeReactionField(system.Mmn.auxsize());

  const votca::Index qptotal = opt.qpmax - opt.qpmin + 1;
  Eigen::MatrixXd ref = Eigen::MatrixXd::Zero(qptotal, qptotal);
  for (votca::Index n1 = 0; n1 < qptotal; ++n1) {
    for (votca::Index n2 = 0; n2 < qptotal; ++n2) {
      for (votca::Index m = 0; m < system.Mmn.nsize(); ++m) {
        const double s = (m + opt.rpamin <= opt.homo) ? -1.0 : 1.0;
        ref(n1, n2) +=
            0.5 * s *
            double(system.Mmn[n1 + opt.qpmin - opt.rpamin].row(m) * R *
                   system.Mmn[n2 + opt.qpmin - opt.rpamin].row(m).transpose());
      }
    }
  }

  GW gw(system.log, system.Mmn, system.vxc, system.mo_eigenvalues);
  gw.configure(opt);
  gw.setReactionField(R);
  gw.CalculateGWPerturbation();
  const Eigen::MatrixXd& sigma = gw.getSigmaReac();
  BOOST_REQUIRE_EQUAL(sigma.rows(), qptotal);
  BOOST_CHECK_SMALL((sigma - ref).cwiseAbs().maxCoeff(), 1e-12);
  BOOST_CHECK_SMALL((sigma - sigma.transpose()).cwiseAbs().maxCoeff(), 1e-14);
  libint2::finalize();
}

// The PPM rotates the auxiliary index of Mmn in place. The bare Coulomb
// terms do not notice; a reaction field does, and must be brought along.
// A second GW run on the same, now rotated, integrals has to give the same
// reaction self-energy and the same QP energies as the first.
BOOST_AUTO_TEST_CASE(reaction_field_follows_rotated_aux_frame) {
  if (!libint2::initialized()) libint2::initialize();
  GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
  const GW::options opt = MakeGWTestOptions();
  const Eigen::MatrixXd R = SomeReactionField(system.Mmn.auxsize());

  GW gw1(system.log, system.Mmn, system.vxc, system.mo_eigenvalues);
  gw1.configure(opt);
  gw1.setReactionField(R);
  gw1.CalculateGWPerturbation();
  BOOST_REQUIRE_GT(system.Mmn.AuxFrame().size(), 0);  // it was rotated

  GW gw2(system.log, system.Mmn, system.vxc, system.mo_eigenvalues);
  gw2.configure(opt);
  gw2.setReactionField(R);
  gw2.CalculateGWPerturbation();
  BOOST_CHECK_SMALL(
      (gw1.getSigmaReac() - gw2.getSigmaReac()).cwiseAbs().maxCoeff(), 1e-12);
  BOOST_CHECK_SMALL(
      (gw1.getGWAResults() - gw2.getGWAResults()).cwiseAbs().maxCoeff(), 1e-8);

  // And it is actually there: without it the QP energies differ.
  GW gw0(system.log, system.Mmn, system.vxc, system.mo_eigenvalues);
  gw0.configure(opt);
  gw0.CalculateGWPerturbation();
  BOOST_CHECK_GT(
      (gw0.getGWAResults() - gw1.getGWAResults()).cwiseAbs().maxCoeff(), 1e-4);
  libint2::finalize();
}

// QSGW snapshots Mmn through operator[] and writes it back every iteration.
// Whether the snapshot is taken in the fill-time frame (integrals refilled
// after the seed run) or in the frame the seed's PPM left behind must not
// matter, with a reaction field as without.
BOOST_AUTO_TEST_CASE(qsgw_reaction_field_independent_of_aux_frame) {
  if (!libint2::initialized()) libint2::initialize();
  GW::options opt = MakeGWTestOptions();
  opt.do_qsgw = true;
  opt.qsgw_max_iterations = 50;
  opt.qsgw_sc_limit = 1e-6;
  opt.gw_mixing_order = 20;
  opt.gw_mixing_alpha = 0.2;
  opt.qpmax = 13;
  opt.qsgw_max_virt_correction = 0.35;
  const votca::Index qptotal = opt.qpmax - opt.qpmin + 1;

  auto run = [&](bool refill) {
    GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
    const Eigen::MatrixXd vxc_gw =
        system.vxc.block(opt.qpmin, opt.qpmin, qptotal, qptotal);
    const Eigen::MatrixXd R = SomeReactionField(system.Mmn.auxsize());
    GW gw(system.log, system.Mmn, vxc_gw, system.mo_eigenvalues);
    gw.configure(opt);
    gw.setReactionField(R);
    gw.CalculateGWPerturbation();
    if (refill) {
      // Not Rebuild(): GWTestSystem's basis objects are gone by now.
      BasisSet basis;
      basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
      QMMolecule mol("methane", 0);
      mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/gw/molecule.xyz");
      AOBasis aobasis;
      aobasis.Fill(basis, mol);
      system.Mmn.Fill(aobasis, aobasis, system.mo_eigenvectors);
    }
    gw.CalculateQSGW();
    return Eigen::VectorXd(gw.getGWAResults());
  };
  const Eigen::VectorXd e_refilled = run(true);
  const Eigen::VectorXd e_rotated = run(false);
  BOOST_CHECK_SMALL((e_refilled - e_rotated).cwiseAbs().maxCoeff(), 1e-6);
  libint2::finalize();
}

// Step 5. With R = 0 the dressing is the identity, and everything --
// exchange, RPA, Sigma_c, the QP energies -- must be what it is without an
// environment at all.
BOOST_AUTO_TEST_CASE(zero_reaction_field_changes_nothing) {
  if (!libint2::initialized()) libint2::initialize();
  const GW::options opt = MakeGWTestOptions();
  GWTestSystem plain("mo_eigenvectors.mm", "vxc.mm");
  const Eigen::VectorXd e_plain = RunGWPerturbation(plain, opt);

  GWTestSystem zero("mo_eigenvectors.mm", "vxc.mm");
  GW gw(zero.log, zero.Mmn, zero.vxc, zero.mo_eigenvalues);
  gw.configure(opt);
  gw.setReactionField(
      Eigen::MatrixXd::Zero(zero.Mmn.auxsize(), zero.Mmn.auxsize()));
  gw.CalculateGWPerturbation();
  BOOST_CHECK(zero.Mmn.Dressed());
  BOOST_CHECK_SMALL((gw.getGWAResults() - e_plain).cwiseAbs().maxCoeff(),
                    1e-10);
  libint2::finalize();
}

// Dressing the auxiliary index must not touch the exchange: it is with the
// bare v, in whatever frame the integrals are in -- here after an arbitrary
// orthogonal rotation, as the PPM would leave them, then dressed.
BOOST_AUTO_TEST_CASE(exchange_is_bare_on_dressed_integrals) {
  if (!libint2::initialized()) libint2::initialize();
  GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
  const GW::options gwopt = MakeGWTestOptions();
  RPA rpa(system.log, system.Mmn);
  rpa.configure(gwopt.homo, gwopt.rpamin, gwopt.rpamax);
  std::unique_ptr<Sigma_base> sigma =
      SigmaFactory().Create("ppm", system.Mmn, rpa);
  Sigma_base::options opt;
  opt.homo = gwopt.homo;
  opt.qpmin = gwopt.qpmin;
  opt.qpmax = gwopt.qpmax;
  opt.rpamin = gwopt.rpamin;
  opt.rpamax = gwopt.rpamax;
  opt.eta = gwopt.eta;
  opt.order = 0;
  opt.alpha = 0.0;
  sigma->configure(opt);
  const Eigen::MatrixXd sx_bare = sigma->CalcExchangeMatrix();

  const votca::Index n = system.Mmn.auxsize();
  std::srand(7);
  const Eigen::HouseholderQR<Eigen::MatrixXd> qr(
      Eigen::MatrixXd(Eigen::MatrixXd::Random(n, n)));
  const Eigen::MatrixXd Q = qr.householderQ();
  system.Mmn.MultiplyRightWithAuxMatrix(Q);
  BOOST_CHECK(system.Mmn.AuxFrameIsOrthogonal());
  BOOST_CHECK_SMALL(
      (sigma->CalcExchangeMatrix() - sx_bare).cwiseAbs().maxCoeff(), 1e-12);

  const Eigen::MatrixXd R = SomeReactionField(n);
  system.Mmn.DressAuxIndex(EnvironmentScreening::DressingMatrix(R));
  BOOST_CHECK(!system.Mmn.AuxFrameIsOrthogonal());
  BOOST_CHECK_SMALL(
      (sigma->CalcExchangeMatrix() - sx_bare).cwiseAbs().maxCoeff(), 1e-11);
  libint2::finalize();
}

// With mixing, the evGW energies must solve the QP equation with the
// screening of the last iteration, E = e_DFT + Sigma_x - V_xc + Sigma_c(E),
// not be Sigma_c at the mixed energies. Two iterations: the screening of the
// second is built from the G0W0 energies.
BOOST_AUTO_TEST_CASE(evgw_energies_solve_their_qp_equation) {
  if (!libint2::initialized()) libint2::initialize();
  GW::options opt = MakeGWTestOptions();
  opt.gw_mixing_order = 20;
  opt.gw_mixing_alpha = 0.5;
  opt.g_sc_limit = 1e-8;

  GWTestSystem g0w0_system("mo_eigenvectors.mm", "vxc.mm");
  const Eigen::VectorXd e_g0w0 = RunGWPerturbation(g0w0_system, opt);
  opt.gw_sc_max_iterations = 2;
  GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
  const Eigen::VectorXd E = RunGWPerturbation(system, opt);

  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  QMMolecule mol("methane", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/gw/molecule.xyz");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  TCMatrix_gwbse M;
  M.Initialize(aobasis.AOBasisSize(), opt.rpamin, opt.rpamax, opt.rpamin,
               opt.rpamax);
  M.Fill(aobasis, aobasis, system.mo_eigenvectors);
  Logger log;
  RPA rpa(log, M);
  rpa.configure(opt.homo, opt.rpamin, opt.rpamax);
  rpa.UpdateRPAInputEnergies(system.mo_eigenvalues, e_g0w0, opt.qpmin);
  std::unique_ptr<Sigma_base> sigma =
      SigmaFactory().Create(opt.sigma_integration, M, rpa);
  Sigma_base::options sopt;
  sopt.homo = opt.homo;
  sopt.qpmin = opt.qpmin;
  sopt.qpmax = opt.qpmax;
  sopt.rpamin = opt.rpamin;
  sopt.rpamax = opt.rpamax;
  sopt.eta = opt.eta;
  sigma->configure(sopt);
  const Eigen::VectorXd sigma_x = sigma->CalcExchangeMatrix().diagonal();
  sigma->PrepareScreening();
  const Eigen::VectorXd residual =
      system.mo_eigenvalues.head(E.size()) + sigma_x -
      system.vxc.diagonal().head(E.size()) + sigma->CalcCorrelationDiag(E) - E;
  BOOST_TEST_MESSAGE("evGW QP equation residual "
                     << residual.cwiseAbs().maxCoeff());
  BOOST_CHECK_SMALL(residual.cwiseAbs().maxCoeff(), 1e-6);
  libint2::finalize();
}

// Adaptive screening updates (G iterated with W fixed, W rebuilt when the
// energies have settled for it) reach the same evGW solution as rebuilding
// W in every iteration, with fewer screening builds: the energies solve the
// QP equation with W built from themselves.
BOOST_AUTO_TEST_CASE(evgw_adaptive_screening_update_same_fixed_point) {
  if (!libint2::initialized()) libint2::initialize();
  for (const std::string integration : {"ppm", "exact"}) {
    GW::options opt = MakeGWTestOptions();
    opt.sigma_integration = integration;
    opt.reset_3c = 0;
    opt.gw_sc_max_iterations = 100;
    opt.gw_sc_limit = 1e-7;
    opt.g_sc_limit = 1e-10;
    // a short history: with 5 or 20 the classic scheme stalls at 1e-5..1e-6
    // on the top level with exact Sigma here (the adaptive one does not)
    opt.gw_mixing_order = 3;
    opt.gw_mixing_alpha = 0.7;

    auto run = [&](const std::string& mode, votca::Index& builds,
                   votca::Index& iters) {
      GW::options o = opt;
      o.screening_update = mode;
      GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
      GW gw(system.log, system.Mmn, system.vxc, system.mo_eigenvalues);
      gw.configure(o);
      gw.CalculateGWPerturbation();
      builds = gw.ScreeningBuilds();
      iters = gw.QPIterations();
      return Eigen::VectorXd(gw.getGWAResults());
    };
    votca::Index builds_every = 0, iters_every = 0, builds_adaptive = 0,
                 iters_adaptive = 0;
    const Eigen::VectorXd e_every = run("every", builds_every, iters_every);
    const Eigen::VectorXd e_adaptive =
        run("adaptive", builds_adaptive, iters_adaptive);
    BOOST_TEST_MESSAGE(integration << ": every " << iters_every
                                   << " iterations / " << builds_every
                                   << " builds, adaptive " << iters_adaptive
                                   << " / " << builds_adaptive);
    BOOST_CHECK_EQUAL(builds_every, iters_every);
    BOOST_CHECK_LT(builds_adaptive, builds_every);
    BOOST_CHECK_SMALL((e_every - e_adaptive).cwiseAbs().maxCoeff(), 1e-6);

    // the adaptive energies solve the QP equation with W from themselves
    GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
    Logger log;
    RPA rpa(log, system.Mmn);
    rpa.configure(opt.homo, opt.rpamin, opt.rpamax);
    rpa.UpdateRPAInputEnergies(system.mo_eigenvalues, e_adaptive, opt.qpmin);
    std::unique_ptr<Sigma_base> sigma =
        SigmaFactory().Create(opt.sigma_integration, system.Mmn, rpa);
    Sigma_base::options sopt;
    sopt.homo = opt.homo;
    sopt.qpmin = opt.qpmin;
    sopt.qpmax = opt.qpmax;
    sopt.rpamin = opt.rpamin;
    sopt.rpamax = opt.rpamax;
    sopt.eta = opt.eta;
    sigma->configure(sopt);
    const Eigen::VectorXd sigma_x = sigma->CalcExchangeMatrix().diagonal();
    sigma->PrepareScreening();
    const Eigen::VectorXd residual =
        system.mo_eigenvalues.head(e_adaptive.size()) + sigma_x -
        system.vxc.diagonal().head(e_adaptive.size()) +
        sigma->CalcCorrelationDiag(e_adaptive) - e_adaptive;
    BOOST_CHECK_SMALL(residual.cwiseAbs().maxCoeff(), 1e-6);
  }
  libint2::finalize();
}

// Root tracking (local QP search near the previous root after the first
// iteration, verified by a full search after convergence) gives the same
// evGW energies as the full search in every iteration, with and without
// adaptive screening updates, and actually solves levels locally.
BOOST_AUTO_TEST_CASE(evgw_root_tracking_same_result) {
  if (!libint2::initialized()) libint2::initialize();
  for (const std::string integration : {"ppm", "exact"}) {
    for (const std::string update : {"every", "adaptive"}) {
      GW::options opt = MakeGWTestOptions();
      opt.sigma_integration = integration;
      opt.screening_update = update;
      opt.reset_3c = 0;
      opt.gw_sc_max_iterations = 100;
      opt.gw_sc_limit = 1e-7;
      opt.g_sc_limit = 1e-10;
      opt.gw_mixing_order = 3;
      opt.gw_mixing_alpha = 0.7;
      auto run = [&](bool tracking, votca::Index& local, votca::Index& iters) {
        GW::options o = opt;
        o.qp_root_tracking = tracking;
        GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
        GW gw(system.log, system.Mmn, system.vxc, system.mo_eigenvalues);
        gw.configure(o);
        gw.CalculateGWPerturbation();
        local = gw.TrackedLocalSolves();
        iters = gw.QPIterations();
        return Eigen::VectorXd(gw.getGWAResults());
      };
      votca::Index local_off = 0, iters_off = 0, local_on = 0, iters_on = 0;
      const Eigen::VectorXd e_full = run(false, local_off, iters_off);
      const Eigen::VectorXd e_tracked = run(true, local_on, iters_on);
      BOOST_TEST_MESSAGE(integration << "/" << update << ": full search "
                                     << iters_off << " iterations, tracking "
                                     << iters_on << " iterations, " << local_on
                                     << " local level solves");
      BOOST_CHECK_EQUAL(local_off, 0);
      BOOST_CHECK_GT(local_on, 0);
      BOOST_CHECK_SMALL((e_full - e_tracked).cwiseAbs().maxCoeff(), 1e-6);
    }
  }
  libint2::finalize();
}

// GW_UKS on a closed shell (alpha = beta) is evGW of GW: same energies for
// both spins as the restricted code, also with adaptive screening updates
// and root tracking, which in addition must reproduce the "every" / full
// search energies of GW_UKS itself.
namespace {
GW_UKS::options ToUKS(const GW::options& opt) {
  GW_UKS::options u;
  u.homo_alpha = opt.homo;
  u.homo_beta = opt.homo;
  u.qpmin = opt.qpmin;
  u.qpmax = opt.qpmax;
  u.rpamin = opt.rpamin;
  u.rpamax = opt.rpamax;
  u.eta = opt.eta;
  u.g_sc_limit = opt.g_sc_limit;
  u.g_sc_max_iterations = opt.g_sc_max_iterations;
  u.gw_sc_limit = opt.gw_sc_limit;
  u.gw_sc_max_iterations = opt.gw_sc_max_iterations;
  u.shift = opt.shift;
  u.ScaHFX = opt.ScaHFX;
  u.sigma_integration = opt.sigma_integration;
  u.reset_3c = opt.reset_3c;
  u.qp_solver = opt.qp_solver;
  u.qp_solver_alpha = opt.qp_solver_alpha;
  u.qp_grid_steps = opt.qp_grid_steps;
  u.qp_grid_spacing = opt.qp_grid_spacing;
  u.qp_full_window_half_width = opt.qp_full_window_half_width;
  u.qp_dense_spacing = opt.qp_dense_spacing;
  u.qp_adaptive_shell_width = opt.qp_adaptive_shell_width;
  u.qp_adaptive_shell_count = opt.qp_adaptive_shell_count;
  u.gw_mixing_order = opt.gw_mixing_order;
  u.gw_mixing_alpha = opt.gw_mixing_alpha;
  u.quadrature_scheme = opt.quadrature_scheme;
  u.order = opt.order;
  u.alpha = opt.alpha;
  u.out_of_window_shift = opt.out_of_window_shift;
  u.qp_restrict_search = opt.qp_restrict_search;
  u.qp_root_continuity = opt.qp_root_continuity;
  u.screening_update = opt.screening_update;
  u.qp_root_tracking = opt.qp_root_tracking;
  u.screening_update_ratio = opt.screening_update_ratio;
  u.screening_update_max_inner = opt.screening_update_max_inner;
  u.qp_zero_margin = opt.qp_zero_margin;
  u.qp_virtual_min_energy = opt.qp_virtual_min_energy;
  u.qp_root_finder = opt.qp_root_finder;
  u.qp_grid_search_mode = opt.qp_grid_search_mode;
  return u;
}

struct UKSRun {
  Eigen::VectorXd alpha;
  Eigen::VectorXd beta;
  votca::Index builds = 0;
  votca::Index iterations = 0;
  votca::Index local = 0;
};

UKSRun RunEvGW_UKS(const GW::options& opt) {
  GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
  TCMatrix_gwbse_spin Mmn;
  Mmn.alpha = system.Mmn;
  Mmn.beta = system.Mmn;
  GW_UKS gw(system.log, Mmn, system.vxc, system.vxc, system.mo_eigenvalues,
            system.mo_eigenvalues);
  gw.configure(ToUKS(opt));
  gw.CalculateGWPerturbation();
  UKSRun r;
  r.alpha = gw.getGWAResultsAlpha();
  r.beta = gw.getGWAResultsBeta();
  r.builds = gw.ScreeningBuilds();
  r.iterations = gw.QPIterations();
  r.local = gw.TrackedLocalSolves();
  return r;
}
}  // namespace

BOOST_AUTO_TEST_CASE(evgw_uks_closed_shell_equals_rks) {
  if (!libint2::initialized()) libint2::initialize();
  for (const std::string integration : {"ppm", "exact"}) {
    GW::options opt = MakeGWTestOptions();
    opt.sigma_integration = integration;
    opt.reset_3c = 0;
    opt.gw_sc_max_iterations = 100;
    opt.gw_sc_limit = 1e-7;
    opt.g_sc_limit = 1e-10;
    opt.gw_mixing_order = 3;
    opt.gw_mixing_alpha = 0.7;

    const Eigen::VectorXd e_rks = [&]() {
      GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
      return RunGWPerturbation(system, opt);
    }();

    const UKSRun every = RunEvGW_UKS(opt);
    BOOST_CHECK_EQUAL(every.builds, every.iterations);
    BOOST_CHECK_SMALL((every.alpha - every.beta).cwiseAbs().maxCoeff(), 1e-10);
    BOOST_CHECK_SMALL((every.alpha - e_rks).cwiseAbs().maxCoeff(), 1e-6);

    GW::options o = opt;
    o.screening_update = "adaptive";
    const UKSRun adaptive = RunEvGW_UKS(o);
    BOOST_CHECK_LT(adaptive.builds, every.builds);
    BOOST_CHECK_EQUAL(adaptive.local, 0);
    BOOST_CHECK_SMALL((adaptive.alpha - e_rks).cwiseAbs().maxCoeff(), 1e-6);
    BOOST_CHECK_SMALL((adaptive.beta - e_rks).cwiseAbs().maxCoeff(), 1e-6);

    for (const std::string update : {"every", "adaptive"}) {
      GW::options t = opt;
      t.screening_update = update;
      t.qp_root_tracking = true;
      const UKSRun tracked = RunEvGW_UKS(t);
      BOOST_TEST_MESSAGE("UKS " << integration << "/" << update << " tracking: "
                                << tracked.iterations << " iterations, "
                                << tracked.builds << " builds, "
                                << tracked.local << " local level solves");
      BOOST_CHECK_GT(tracked.local, 0);
      BOOST_CHECK_SMALL((tracked.alpha - e_rks).cwiseAbs().maxCoeff(), 1e-6);
      BOOST_CHECK_SMALL((tracked.beta - e_rks).cwiseAbs().maxCoeff(), 1e-6);
    }
    BOOST_TEST_MESSAGE("UKS " << integration << ": every " << every.iterations
                              << " iterations, adaptive " << adaptive.iterations
                              << " / " << adaptive.builds << " builds");
  }
  libint2::finalize();
}

// Batched UKS Sigma_c (as used by the QP scan through Prefetch) equals the
// element-wise sums: in chunks too small for the Chebyshev split to 1e-13,
// and single frequencies to 1e-11 (rounded differently near poles, see
// test_sigma_exact.cc).
BOOST_AUTO_TEST_CASE(sigma_uks_batched_equals_direct) {
  if (!libint2::initialized()) libint2::initialize();
  for (const std::string integration : {"ppm", "exact"}) {
    GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
    TCMatrix_gwbse_spin Mmn;
    Mmn.alpha = system.Mmn;
    Mmn.beta = system.Mmn;
    Logger log;
    RPA_UKS rpa(log, Mmn);
    rpa.configure(4, 4, 0, 16);
    rpa.setRPAInputEnergies(system.mo_eigenvalues, system.mo_eigenvalues);
    PPM ppm;
    std::unique_ptr<Sigma_base_UKS> sigma = SigmaFactory_UKS().Create(
        integration, Mmn, rpa, TCMatrix::SpinChannel::Beta);
    if (integration == "ppm") {
      ppm.PPM_construct_parameters(rpa);
      dynamic_cast<Sigma_PPM_UKS&>(*sigma).SetSharedPPM(ppm);
    }
    Sigma_base_UKS::options opt;
    opt.homo = 4;
    opt.qpmin = 0;
    opt.qpmax = 16;
    opt.rpamin = 0;
    opt.rpamax = 16;
    opt.eta = 1e-3;
    sigma->configure(opt);
    sigma->PrepareScreening();
    for (votca::Index level = 0; level < 17; ++level) {
      for (double halfwidth : {0.75, 3.0}) {
        const votca::Index n = (halfwidth < 1.0) ? 63 : 1001;
        const double e = system.mo_eigenvalues(level);
        const Eigen::VectorXd w =
            Eigen::VectorXd::LinSpaced(n, e - halfwidth, e + halfwidth);
        const Eigen::VectorXd batched =
            sigma->CalcCorrelationDiagElements(level, w);
        double cheb = 0.0;
        double single = 0.0;
        double scale = 0.0;
        for (votca::Index i0 = 0; i0 < n; i0 += 16) {
          const votca::Index len = std::min<votca::Index>(16, n - i0);
          const Eigen::VectorXd chunk =
              sigma->CalcCorrelationDiagElements(level, w.segment(i0, len));
          for (votca::Index i = 0; i < len; ++i) {
            const double direct =
                sigma->CalcCorrelationDiagElement(level, w(i0 + i));
            cheb = std::max(cheb, std::abs(batched(i0 + i) - chunk(i)));
            single = std::max(single, std::abs(batched(i0 + i) - direct));
            scale = std::max(scale, std::abs(direct));
          }
        }
        BOOST_CHECK_SMALL(cheb, 1e-13 * std::max(1.0, scale));
        BOOST_CHECK_SMALL(single, 1e-11 * std::max(1.0, scale));
      }
    }
  }
  libint2::finalize();
}

// Frozen-core Sigma_x (Mmn from rpamin = 1 on) plus the exchange with the
// frozen core level equals the all-electron Sigma_x in the window.
BOOST_AUTO_TEST_CASE(core_exchange_completes_frozen_core_sigma_x) {
  if (!libint2::initialized()) libint2::initialize();
  GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/gw/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  auto sigma_x = [&](TCMatrix_gwbse& Mmn, votca::Index rpamin) {
    Logger log;
    RPA rpa(log, Mmn);
    rpa.configure(4, rpamin, 16);
    rpa.setRPAInputEnergies(system.mo_eigenvalues.segment(rpamin, 17 - rpamin));
    std::unique_ptr<Sigma_base> sigma = SigmaFactory().Create("ppm", Mmn, rpa);
    Sigma_base::options opt;
    opt.homo = 4;
    opt.qpmin = 1;
    opt.qpmax = 16;
    opt.rpamin = rpamin;
    opt.rpamax = 16;
    opt.eta = 1e-3;
    sigma->configure(opt);
    return Eigen::MatrixXd(sigma->CalcExchangeMatrix());
  };
  // all electron: the system's Mmn (levels 0..16)
  const Eigen::MatrixXd x_all = sigma_x(system.Mmn, 0);
  TCMatrix_gwbse Mmn_fc;
  Mmn_fc.Initialize(aobasis.AOBasisSize(), 1, 16, 1, 16);
  Mmn_fc.Fill(aobasis, aobasis, system.mo_eigenvectors);
  const Eigen::MatrixXd x_fc = sigma_x(Mmn_fc, 1);
  const Eigen::MatrixXd core = TCMatrix_gwbse::CoreExchange(
      aobasis, aobasis, system.mo_eigenvectors, 1, 1, 16);
  BOOST_CHECK_GT((x_all - x_fc).cwiseAbs().maxCoeff(), 1e-3);
  BOOST_CHECK_SMALL((x_fc + core - x_all).cwiseAbs().maxCoeff(), 1e-10);
  libint2::finalize();
}

// The QSGW fixed point, checked from scratch in the converged QP basis
// phi = psi U: with the integrals refilled for those orbitals and the QP
// energies E in the RPA, H = U^T H0 U + Sigma_x + 1/2 [Sigma_c(E_i) +
// Sigma_c(E_j)] must be diagonal with eigenvalues E. qpmin > rpamin, so the
// core level is in the screening but not in the QP window.
namespace {
double QSGWFixedPointResidual(const std::string& sigma_integration) {
  GW::options opt = MakeGWTestOptions();
  opt.sigma_integration = sigma_integration;
  opt.do_qsgw = true;
  opt.qsgw_max_iterations = 40;
  opt.qsgw_sc_limit = 1e-9;
  opt.gw_mixing_order = 20;
  opt.gw_mixing_alpha = 0.7;
  opt.qpmin = 1;
  opt.qpmax = 13;
  opt.qsgw_max_virt_correction = 10.0;  // no trimming
  const votca::Index qptotal = opt.qpmax - opt.qpmin + 1;

  GWTestSystem system("mo_eigenvectors.mm", "vxc.mm");
  const Eigen::MatrixXd vxc_gw =
      system.vxc.block(opt.qpmin, opt.qpmin, qptotal, qptotal);
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  QMMolecule mol("methane", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/gw/molecule.xyz");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);

  GW gw(system.log, system.Mmn, vxc_gw, system.mo_eigenvalues);
  gw.configure(opt);
  gw.CalculateGWPerturbation();
  system.Mmn.Fill(aobasis, aobasis, system.mo_eigenvectors);
  gw.CalculateQSGW();
  const Eigen::VectorXd E = gw.getGWAResults();
  const Eigen::MatrixXd U = gw.getQSGWRotation();
  {
    std::stringstream log_text;
    log_text << system.log;
    const std::string text = log_text.str();
    const std::size_t pos = text.find("QSGW converged in");
    BOOST_REQUIRE(pos != std::string::npos);
    BOOST_TEST_MESSAGE(sigma_integration
                       << ": " << text.substr(pos, text.find('\n', pos) - pos));
  }

  Eigen::MatrixXd C_qp = system.mo_eigenvectors;
  C_qp.middleCols(opt.qpmin, qptotal) =
      system.mo_eigenvectors.middleCols(opt.qpmin, qptotal) * U;
  TCMatrix_gwbse M;
  M.Initialize(aobasis.AOBasisSize(), opt.rpamin, opt.rpamax, opt.rpamin,
               opt.rpamax);
  M.Fill(aobasis, aobasis, C_qp);
  Logger log;
  RPA rpa(log, M);
  rpa.configure(opt.homo, opt.rpamin, opt.rpamax);
  rpa.UpdateRPAInputEnergies(system.mo_eigenvalues, E, opt.qpmin);
  std::unique_ptr<Sigma_base> sigma =
      SigmaFactory().Create(opt.sigma_integration, M, rpa);
  Sigma_base::options sopt;
  sopt.homo = opt.homo;
  sopt.qpmin = opt.qpmin;
  sopt.qpmax = opt.qpmax;
  sopt.rpamin = opt.rpamin;
  sopt.rpamax = opt.rpamax;
  sopt.eta = opt.eta;
  sigma->configure(sopt);
  sigma->PrepareScreening();
  Eigen::MatrixXd H0 = -vxc_gw;
  H0.diagonal() += system.mo_eigenvalues.segment(opt.qpmin, qptotal);
  Eigen::MatrixXd H = U.transpose() * H0 * U + sigma->CalcExchangeMatrix() +
                      sigma->CalcCorrelationOffDiag(E);
  H.diagonal() += sigma->CalcCorrelationDiag(E);
  H.diagonal() -= E;
  return H.cwiseAbs().maxCoeff();
}
}  // namespace

BOOST_AUTO_TEST_CASE(qsgw_fixed_point_in_qp_basis) {
  if (!libint2::initialized()) libint2::initialize();
  for (const std::string sigma : {"ppm", "exact"}) {
    const double residual = QSGWFixedPointResidual(sigma);
    BOOST_TEST_MESSAGE("QSGW fixed-point residual (" << sigma
                                                     << "): " << residual);
    BOOST_CHECK_SMALL(residual, 1e-6);
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
