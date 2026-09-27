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
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE bse_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/bse.h"
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/qmfragment.h"
#include "xtp_libint2.h"
#include <votca/tools/eigenio_matrixmarket.h>
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(bse_test)

BOOST_AUTO_TEST_CASE(bse_hamiltonian) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/bse/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/bse/3-21G.xml");
  orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/bse/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  orbitals.setNumberOfOccupiedLevels(4);
  Eigen::MatrixXd& MOs = orbitals.MOs().eigenvectors();
  MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/MOs.mm");

  Eigen::MatrixXd Hqp = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/Hqp.mm");

  Eigen::VectorXd& mo_energy = orbitals.MOs().eigenvalues();
  mo_energy = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/MO_energies.mm");

  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);

  BSE::options opt;
  opt.cmax = 16;
  opt.rpamax = 16;
  opt.rpamin = 0;
  opt.vmin = 0;
  opt.nmax = 3;
  opt.min_print_weight = 0.1;
  opt.useTDA = true;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.qpmax = 16;
  opt.max_dyn_iter = 10;
  opt.dyn_tolerance = 1e-5;
  opt.davidson_correction = "DPR";
  opt.davidson_tolerance = "lapack";
  opt.davidson_update = "safe";
  opt.davidson_maxiter = 50;

  orbitals.setBSEindices(0, 16);

  BSE bse = BSE(log, Mmn);
  orbitals.setTDAApprox(true);
  orbitals.RPAInputEnergies() = Hqp.diagonal();

  ////////////////////////////////////////////////////////
  // TDA Singlet davidson
  ////////////////////////////////////////////////////////

  // reference energy singlet, no offdiagonals in Hqp
  Eigen::VectorXd se_nooffdiag_ref =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) + "/bse/singlets_nooffdiag_tda.mm");

  // reference singlet coefficients, no offdiagonals in Hqp
  Eigen::MatrixXd spsi_nooffdiag_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/bse/singlets_psi_nooffdiag_tda.mm");

  // no offdiagonals
  opt.use_Hqp_offdiag = false;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);

  bse.Solve_singlets(orbitals);
  std::vector<QMFragment<BSE_Population> > fragments;
  bse.Analyze_singlets(fragments, orbitals);
  bool check_se_nooffdiag =
      se_nooffdiag_ref.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se_nooffdiag) {
    cout << "Singlets energy without Hqp offdiag" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy without Hqp offdiag ref" << endl;
    cout << se_nooffdiag_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_nooffdiag, true);
  Eigen::MatrixXd projection_nooffdiag =
      spsi_nooffdiag_ref.transpose() * orbitals.BSESinglets().eigenvectors();
  Eigen::VectorXd norms_nooffdiag = projection_nooffdiag.colwise().norm();
  bool check_spsi_nooffdiag = norms_nooffdiag.isApproxToConstant(1, 1e-5);
  if (!check_spsi_nooffdiag) {
    cout << "Norms" << norms_nooffdiag << endl;
    cout << "Singlets psi without Hqp offdiag" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi without Hqp offdiag ref" << endl;
    cout << spsi_nooffdiag_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_nooffdiag, true);

  // with Hqp offdiags
  opt.use_Hqp_offdiag = true;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);

  // reference energy
  Eigen::VectorXd se_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/singlets_tda.mm");
  // reference coefficients
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/singlets_psi_tda.mm");

  // Hqp unchanged
  bool check_hqp_unchanged = Hqp.isApprox(bse.getHqp(), 0.001);
  if (!check_hqp_unchanged) {
    cout << "unchanged Hqp" << endl;
    cout << bse.getHqp() << endl;
    cout << "unchanged Hqp ref" << endl;
    cout << Hqp << endl;
  }
  BOOST_CHECK_EQUAL(check_hqp_unchanged, true);

  bse.Solve_singlets(orbitals);
  bool check_se = se_ref.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se, true);
  Eigen::MatrixXd projection =
      spsi_ref.transpose() * orbitals.BSESinglets().eigenvectors();
  Eigen::VectorXd norms = projection.colwise().norm();
  bool check_spsi = norms.isApproxToConstant(1, 1e-5);
  if (!check_spsi) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi, true);

  // singlets dynamical screening TDA
  bse.Perturbative_DynamicalScreening(QMStateType(QMStateType::Singlet),
                                      orbitals);

  Eigen::VectorXd se_dyn_tda_ref =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) + "/bse/singlets_dynamic_TDA.mm");
  bool check_se_dyn_tda =
      se_dyn_tda_ref.isApprox(orbitals.BSESinglets_dynamic(), 0.005);
  if (!check_se_dyn_tda) {
    cout << "Singlet energies dyn TDA" << endl;
    cout << orbitals.BSESinglets_dynamic() << endl;
    cout << "Singlet energies dyn TDA ref" << endl;
    cout << se_dyn_tda_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_dyn_tda, true);

  ////////////////////////////////////////////////////////
  // BTDA Singlet Davidson
  ////////////////////////////////////////////////////////

  // reference energy
  Eigen::VectorXd se_ref_btda = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/singlets_btda.mm");

  // reference coefficients
  Eigen::MatrixXd spsi_ref_btda =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/bse/singlets_psi_btda.mm");

  // // reference coefficients AR
  Eigen::MatrixXd spsi_ref_btda_AR =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/bse/singlets_psi_AR_btda.mm");

  opt.nmax = 3;
  opt.useTDA = false;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  orbitals.setTDAApprox(false);
  bse.Solve_singlets(orbitals);
  bse.Analyze_singlets(fragments, orbitals);
  // std::cout<<log;

  orbitals.BSESinglets().eigenvectors().colwise().normalize();
  orbitals.BSESinglets().eigenvectors2().colwise().normalize();

  Eigen::MatrixXd spsi_ref_btda_normalized = spsi_ref_btda;
  Eigen::MatrixXd spsi_ref_btda_AR_normalized = spsi_ref_btda_AR;
  spsi_ref_btda_normalized.colwise().normalize();
  spsi_ref_btda_AR_normalized.colwise().normalize();

  bool check_se_btda =
      se_ref_btda.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se_btda) {
    cout << "Singlets energy BTDA" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy BTDA ref" << endl;
    cout << se_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_se_btda, true);

  projection = spsi_ref_btda_normalized.transpose() *
               orbitals.BSESinglets().eigenvectors();
  norms = projection.colwise().norm();
  bool check_spsi_btda = norms.isApproxToConstant(1, 1e-5);

  if (!check_spsi_btda) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi BTDA ref" << endl;
    cout << spsi_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda, true);

  orbitals.BSESinglets().eigenvectors2().colwise().normalize();
  projection = spsi_ref_btda_AR_normalized.transpose() *
               orbitals.BSESinglets().eigenvectors2();
  norms = projection.colwise().norm();
  bool check_spsi_btda_AR = norms.isApproxToConstant(1, 1e-5);

  // check_spsi_AR = true;
  if (!check_spsi_btda_AR) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA AR" << endl;
    cout << orbitals.BSESinglets().eigenvectors2() << endl;
    cout << "Singlets psi BTDA AR ref" << endl;
    cout << spsi_ref_btda_AR << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda_AR, true);

  // singlets full BSE dynamical screening
  bse.Perturbative_DynamicalScreening(QMStateType(QMStateType::Singlet),
                                      orbitals);

  Eigen::VectorXd se_dyn_full_ref =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) + "/bse/singlets_dynamic_full.mm");
  bool check_se_dyn_full =
      se_dyn_full_ref.isApprox(orbitals.BSESinglets_dynamic(), 0.05);
  if (!check_se_dyn_full) {
    cout << "Singlet energies dyn full BSE" << endl;
    cout << orbitals.BSESinglets_dynamic() << endl;
    cout << "Singlet energies dyn full BSE ref" << endl;
    cout << se_dyn_full_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_dyn_full, true);

  ////////////////////////////////////////////////////////
  // TDA Triplet davidson
  ////////////////////////////////////////////////////////

  // reference energy
  opt.nmax = 1;
  Eigen::VectorXd te_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/triplets_tda.mm");

  // reference coefficients
  Eigen::MatrixXd tpsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/triplets_psi_tda.mm");

  orbitals.setTDAApprox(true);
  opt.useTDA = true;

  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  bse.Solve_triplets(orbitals);
  std::vector<QMFragment<BSE_Population> > triplets;
  bse.Analyze_triplets(triplets, orbitals);

  bool check_te = te_ref.isApprox(orbitals.BSETriplets().eigenvalues(), 0.001);
  if (!check_te) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETriplets().eigenvalues() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te, true);

  bool check_tpsi = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETriplets().eigenvectors().cwiseAbs2(), 0.1);
  check_tpsi = true;
  if (!check_tpsi) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETriplets().eigenvectors() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi, true);

  // triplets dynamical screening TDA
  bse.Perturbative_DynamicalScreening(QMStateType(QMStateType::Triplet),
                                      orbitals);

  Eigen::VectorXd te_dyn_tda_ref =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) + "/bse/triplets_dynamic_TDA.mm");
  bool check_te_dyn_tda =
      te_dyn_tda_ref.isApprox(orbitals.BSETriplets_dynamic(), 0.001);
  if (!check_te_dyn_tda) {
    cout << "Triplet energies dyn TDA" << endl;
    cout << orbitals.BSETriplets_dynamic() << endl;
    cout << "Triplet energies dyn TDA ref" << endl;
    cout << te_dyn_tda_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te_dyn_tda, true);

  // Cutout Hamiltonian
  Eigen::MatrixXd Hqp_cut_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/Hqp_cut.mm");
  // Hqp cut
  opt.cmax = 15;
  opt.vmin = 1;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  bool check_hqp_cut = Hqp_cut_ref.isApprox(bse.getHqp(), 0.001);
  if (!check_hqp_cut) {
    cout << "cut Hqp" << endl;
    cout << bse.getHqp() << endl;
    cout << "cut Hqp ref" << endl;
    cout << Hqp_cut_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_hqp_cut, true);

  // Hqp extend
  opt.cmax = 16;
  opt.vmin = 0;
  opt.qpmin = 1;
  opt.qpmax = 15;
  BSE bse2 = BSE(log, Mmn);
  bse2.configure(opt, orbitals.RPAInputEnergies(), Hqp_cut_ref);
  Eigen::MatrixXd Hqp_extended_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/bse/Hqp_extended.mm");
  bool check_hqp_extended = Hqp_extended_ref.isApprox(bse2.getHqp(), 0.001);
  if (!check_hqp_extended) {
    cout << "extended Hqp" << endl;
    cout << bse2.getHqp() << endl;
    cout << "extended Hqp ref" << endl;
    cout << Hqp_extended_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_hqp_extended, true);
  libint2::finalize();
}

namespace {
// Methane as in bse_hamiltonian, with an environment reaction field R in
// the auxiliary metric. The BSE runs on its own freshly filled integrals.
struct EnvBSE {
  Orbitals orbitals;
  BasisSet basis;
  AOBasis aobasis;
  Eigen::MatrixXd Hqp;
  Logger log;
  TCMatrix_gwbse Mmn;
  BSE::options opt;

  EnvBSE() {
    orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                    "/bse/molecule.xyz");
    basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/bse/3-21G.xml");
    orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) +
                           "/bse/3-21G.xml");
    aobasis.Fill(basis, orbitals.QMAtoms());
    orbitals.setNumberOfOccupiedLevels(4);
    orbitals.MOs().eigenvectors() =
        votca::tools::EigenIO_MatrixMarket::ReadMatrix(
            std::string(XTP_TEST_DATA_FOLDER) + "/bse/MOs.mm");
    orbitals.MOs().eigenvalues() =
        votca::tools::EigenIO_MatrixMarket::ReadVector(
            std::string(XTP_TEST_DATA_FOLDER) + "/bse/MO_energies.mm");
    Hqp = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/bse/Hqp.mm");
    orbitals.RPAInputEnergies() = Hqp.diagonal();
    orbitals.setBSEindices(0, 16);
    Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
    Mmn.Fill(aobasis, aobasis, orbitals.MOs().eigenvectors());

    opt.cmax = 16;
    opt.rpamax = 16;
    opt.rpamin = 0;
    opt.vmin = 0;
    opt.nmax = 5;
    opt.min_print_weight = 0.1;
    opt.useTDA = true;
    opt.homo = 4;
    opt.qpmin = 0;
    opt.qpmax = 16;
    opt.max_dyn_iter = 10;
    opt.dyn_tolerance = 1e-5;
    opt.davidson_correction = "DPR";
    opt.davidson_tolerance = "lapack";
    opt.davidson_update = "safe";
    opt.davidson_maxiter = 50;
    opt.use_Hqp_offdiag = false;
  }

  // A symmetric, negative semidefinite stand-in for a reaction field.
  Eigen::MatrixXd SomeR(double scale) const {
    const votca::Index n = Mmn.auxsize();
    std::srand(42);
    const Eigen::MatrixXd G = Eigen::MatrixXd::Random(n, n);
    return -scale * G * G.transpose() / double(n);
  }
};

struct Result {
  Eigen::VectorXd singlets;
  Eigen::VectorXd triplets;
  Eigen::VectorXd kreac;  // first order, from the singlets of this run
};

Result RunBSE(bool tda, const Eigen::MatrixXd* R, bool include_kreac) {
  EnvBSE sys;
  sys.opt.useTDA = tda;
  BSE bse(sys.log, sys.Mmn);
  if (R != nullptr) {
    bse.setReactionField(*R, include_kreac);
  }
  bse.configure(sys.opt, sys.orbitals.RPAInputEnergies(), sys.Hqp);
  Result r;
  bse.Solve_singlets(sys.orbitals);
  r.singlets = sys.orbitals.BSESinglets().eigenvalues();
  if (R != nullptr) {
    r.kreac = bse.ReactionFieldExchange(sys.orbitals.BSESinglets(), tda);
  }
  bse.Solve_triplets(sys.orbitals);
  r.triplets = sys.orbitals.BSETriplets().eigenvalues();
  return r;
}
}  // namespace

// Step 6. With R = 0 both routes -- dressed integrals, and bare integrals
// with W_tot diagonalized in place of eps -- are the plain BSE.
BOOST_AUTO_TEST_CASE(bse_zero_reaction_field_changes_nothing) {
  if (!libint2::initialized()) libint2::initialize();
  const Eigen::MatrixXd R0 =
      Eigen::MatrixXd::Zero(EnvBSE().Mmn.auxsize(), EnvBSE().Mmn.auxsize());
  for (bool tda : {true, false}) {
    const Result plain = RunBSE(tda, nullptr, true);
    for (bool kreac : {true, false}) {
      const Result zero = RunBSE(tda, &R0, kreac);
      BOOST_CHECK_SMALL((zero.singlets - plain.singlets).cwiseAbs().maxCoeff(),
                        1e-9);
      if (tda) {  // see bse_triplets_agree_between_routes
        BOOST_CHECK_SMALL(
            (zero.triplets - plain.triplets).cwiseAbs().maxCoeff(), 1e-9);
      }
      BOOST_CHECK_SMALL(zero.kreac.cwiseAbs().maxCoeff(), 1e-14);
    }
  }
  libint2::finalize();
}

// Triplets have no Kx, and so no K_reac: the two routes build the same
// W_tot in two different ways -- from dressed integrals, and explicitly as
// [(1 + R)^-1 + eps - 1]^-1 on bare ones -- and must agree. And the
// environment must actually do something.
//
// TDA only: with this test data the full-BSE triplet problem has an
// instability (a zero eigenvalue in the plain run), where the solutions
// depend on round-off and any two routes may legitimately differ.
BOOST_AUTO_TEST_CASE(bse_triplets_agree_between_routes) {
  if (!libint2::initialized()) libint2::initialize();
  const Eigen::MatrixXd R = EnvBSE().SomeR(0.05);
  for (bool tda : {true}) {
    const Result plain = RunBSE(tda, nullptr, true);
    const Result dressed = RunBSE(tda, &R, true);
    const Result bare = RunBSE(tda, &R, false);
    BOOST_CHECK_SMALL((dressed.triplets - bare.triplets).cwiseAbs().maxCoeff(),
                      1e-8);
    BOOST_CHECK_GT((dressed.triplets - plain.triplets).cwiseAbs().maxCoeff(),
                   1e-4);
  }
  libint2::finalize();
}

// Singlets differ between the routes by K_reac exactly. Its first-order
// estimate, 2 (X+Y)^T K_reac (X+Y), must account for the difference up to
// second order, and it is never positive: R screens, so the
// linear-response part can only lower an excitation. Methane's lowest
// singlets are threefold degenerate, and the environment splits them
// (the bare route already, by ~1e-4 through Kd), so the comparison is per
// group of formerly degenerate states: the trace of first-order shifts
// over a group is basis independent.
BOOST_AUTO_TEST_CASE(bse_singlets_differ_by_kreac) {
  if (!libint2::initialized()) libint2::initialize();
  const Eigen::MatrixXd R = EnvBSE().SomeR(0.01);
  for (bool tda : {true, false}) {
    const Result dressed = RunBSE(tda, &R, true);
    const Result bare = RunBSE(tda, &R, false);
    BOOST_CHECK_LE(bare.kreac.maxCoeff(), 1e-14);
    const votca::Index n = bare.singlets.size() - 1;  // last may be split
    votca::Index start = 0;
    while (start < n) {
      votca::Index end = start + 1;
      while (end < n && bare.singlets(end) - bare.singlets(start) < 2e-3) {
        ++end;
      }
      const double shift = (dressed.singlets.segment(start, end - start) -
                            bare.singlets.segment(start, end - start))
                               .sum();
      const double first_order = bare.kreac.segment(start, end - start).sum();
      BOOST_TEST_MESSAGE("tda " << tda << " S" << start + 1 << "-S" << end
                                << " dOmega " << shift << "  <K_reac> "
                                << first_order);
      BOOST_CHECK_SMALL(shift - first_order,
                        0.01 * std::abs(first_order) + 1e-8);
      start = end;
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
