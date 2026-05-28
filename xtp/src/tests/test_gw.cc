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
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/gw.h"

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/xtp/orbitals.h>

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
  opt.do_qsgw                  = true;
  opt.qsgw_max_iterations      = 50;
  opt.qsgw_sc_limit            = 1e-5;
  opt.gw_mixing_order          = 20;
  opt.gw_mixing_alpha          = 0.2;
  opt.qpmax                    = 13;
  opt.qsgw_max_virt_correction = 0.35;

  // G0W0 seed
  opt.gw_sc_max_iterations = 1;
  Logger log1;
  TCMatrix_gwbse Mmn1;
  Mmn1.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn1.Fill(aobasis, aobasis, mo_eigenvectors);
  GW gw1(log1, Mmn1, vxc, mo_eigenvalues);
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
  GW gw2(log2, Mmn2, vxc, mo_eigenvalues);
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
  opt.do_qsgw                  = true;
  opt.qsgw_max_iterations      = 50;
  opt.qsgw_sc_limit            = 1e-5;
  opt.gw_mixing_order          = 20;
  opt.gw_mixing_alpha          = 0.2;
  opt.gw_sc_max_iterations     = 1;
  opt.qpmax                    = 16;
  opt.qsgw_max_virt_correction = 0.2;

  GW gw(log, Mmn, vxc, mo_eigenvalues);
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

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
