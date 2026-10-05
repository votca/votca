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

// Standard includes
#include <algorithm>
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/bse_operator.h"
#include "votca/xtp/bse_operator_uks.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "xtp_libint2.h"
#include <votca/tools/eigenio_matrixmarket.h>
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(bse_operator_test)

BOOST_AUTO_TEST_CASE(bse_operator) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/bse/molecule.xyz");
  orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/bse/3-21G.xml");
  AOBasis aobasis = orbitals.getDftBasis();

  orbitals.setNumberOfOccupiedLevels(4);
  Eigen::MatrixXd& MOs = orbitals.MOs().eigenvectors();
  MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/MOs.mm");

  Eigen::MatrixXd Hqp = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/Hqp.mm");

  Eigen::VectorXd& mo_energy = orbitals.MOs().eigenvalues();
  mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << -0.612601, -0.341755, -0.341755, -0.341755, 0.137304, 0.16678,
      0.16678, 0.16678, 0.671592, 0.671592, 0.671592, 0.974255, 1.01205,
      1.01205, 1.01205, 1.64823, 19.4429;
  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);

  Eigen::MatrixXd rpa_op = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/rpa_op.mm");

  Eigen::VectorXd epsilon_inv = Eigen::VectorXd::Zero(aobasis.AOBasisSize());
  Mmn.MultiplyRightWithAuxMatrix(rpa_op);
  epsilon_inv << 0.999807798016267, 0.994206065211371, 0.917916768047073,
      0.902913813951883, 0.902913745974602, 0.902913584797742,
      0.853352878674581, 0.853352727016914, 0.853352541699637, 0.79703468058566,
      0.797034577207669, 0.797034400395582, 0.787701833916331,
      0.518976361745313, 0.518975064844033, 0.518973712898761,
      0.459286057710524;

  BSEOperator_Options opt;
  opt.cmax = 8;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.rpamin = 0;
  opt.vmin = 0;

  orbitals.setBSEindices(0, 16);
  HqpOperator Hqp_op(epsilon_inv, Mmn, Hqp);
  Hqp_op.configure(opt);
  const Eigen::MatrixXd identity =
      Eigen::MatrixXd::Identity(Hqp_op.rows(), Hqp_op.cols());
  Eigen::MatrixXd hqp_mat = Hqp_op * identity;

  Eigen::MatrixXd hqp_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/hqp_ref.mm");

  bool check_hqp = hqp_mat.isApprox(hqp_ref, 0.001);
  BOOST_CHECK_EQUAL(check_hqp, true);
  bool check_hqpdiag = hqp_mat.diagonal().isApprox(Hqp_op.diagonal(), 0.001);
  BOOST_CHECK_EQUAL(check_hqpdiag, true);
  HxOperator Hx(epsilon_inv, Mmn, Hqp);
  Hx.configure(opt);
  Eigen::MatrixXd hx_mat = Hx * identity;
  Eigen::MatrixXd hx_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/hx_ref.mm");

  bool check_hx = hx_mat.isApprox(hx_ref, 0.001);
  BOOST_CHECK_EQUAL(check_hx, true);
  if (!check_hx) {
    cout << "hx ref" << endl;
    cout << hx_ref << endl;
    cout << "hx result" << endl;
    cout << hx_mat << endl;
  }

  bool check_hxdiag = hx_mat.diagonal().isApprox(Hx.diagonal(), 0.001);
  BOOST_CHECK_EQUAL(check_hxdiag, true);
  HdOperator Hd(epsilon_inv, Mmn, Hqp);
  Hd.configure(opt);
  Eigen::MatrixXd hd_mat = Hd * identity;

  Eigen::MatrixXd hd_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/hd_ref.mm");

  bool check_hd = hd_mat.isApprox(hd_ref, 0.001);

  if (!check_hd) {
    cout << "hd ref" << endl;
    cout << hd_ref << endl;
    cout << "hd result" << endl;
    cout << hd_mat << endl;
  }
  BOOST_CHECK_EQUAL(check_hd, true);

  bool check_hddiag = hd_mat.diagonal().isApprox(Hd.diagonal(), 0.001);
  BOOST_CHECK_EQUAL(check_hddiag, true);

  Hd2Operator Hd2(epsilon_inv, Mmn, Hqp);
  Hd2.configure(opt);
  Eigen::MatrixXd hd2_mat = Hd2 * identity;
  Eigen::MatrixXd hd2_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/hd2_ref.mm");

  bool check_hd2 = hd2_mat.isApprox(hd2_ref, 0.001);
  if (!check_hd2) {
    cout << "hd2 ref" << endl;
    cout << hd2_ref << endl;
    cout << "hd2 result" << endl;
    cout << hd2_mat << endl;
  }
  BOOST_CHECK_EQUAL(check_hd2, true);
  bool check_hd2diag = hd2_mat.diagonal().isApprox(Hd2.diagonal(), 0.001);
  BOOST_CHECK_EQUAL(check_hd2diag, true);

  libint2::finalize();
}

// For a closed-shell reference with identical alpha and beta orbitals, the
// spin-unrestricted BSE must reproduce the restricted singlet and triplet
// spectra: TDA eigenvalues of A, and the full-BSE excitation energies from
// (A-B)(A+B).
BOOST_AUTO_TEST_CASE(uks_closed_shell_equals_singlets_and_triplets) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/bse/molecule.xyz");
  orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/bse/3-21G.xml");
  AOBasis aobasis = orbitals.getDftBasis();
  orbitals.setNumberOfOccupiedLevels(4);
  Eigen::MatrixXd& MOs = orbitals.MOs().eigenvectors();
  MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/MOs.mm");
  const Eigen::MatrixXd Hqp = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/Hqp.mm");

  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);
  const Eigen::MatrixXd rpa_op = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/rpa_op.mm");
  Mmn.MultiplyRightWithAuxMatrix(rpa_op);
  // a strongly screened, non-uniform epsilon^-1, so that W and v differ
  Eigen::VectorXd epsilon_inv(aobasis.AOBasisSize());
  for (votca::Index i = 0; i < epsilon_inv.size(); ++i) {
    epsilon_inv(i) = 0.3 + 0.6 * double(i) / double(epsilon_inv.size());
  }

  BSEOperator_Options opt;
  opt.cmax = 8;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.rpamin = 0;
  opt.vmin = 0;
  BSEOperatorUKS_Options uopt;
  uopt.homo_alpha = opt.homo;
  uopt.homo_beta = opt.homo;
  uopt.rpamin = opt.rpamin;
  uopt.qpmin = opt.qpmin;
  uopt.vmin = opt.vmin;
  uopt.cmax = opt.cmax;

  TCMatrix_gwbse_spin Mspin;
  Mspin.alpha = Mmn;
  Mspin.beta = Mmn;

  auto Dense = [](auto& op) {
    return Eigen::MatrixXd(op *
                           Eigen::MatrixXd::Identity(op.rows(), op.cols()));
  };
  auto Sorted = [](Eigen::VectorXd v) {
    std::sort(v.data(), v.data() + v.size());
    return v;
  };
  auto Concat = [](const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    Eigen::VectorXd c(a.size() + b.size());
    c << a, b;
    return c;
  };
  // squared full-BSE excitation energies, eigenvalues of (A-B)(A+B)
  auto FullBSESquared = [](const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> amb(A - B);
    const Eigen::MatrixXd root = amb.operatorSqrt();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(root * (A + B) * root);
    return Eigen::VectorXd(es.eigenvalues());
  };

  SingletOperator_TDA singlet_A(epsilon_inv, Mmn, Hqp);
  singlet_A.configure(opt);
  TripletOperator_TDA triplet_A(epsilon_inv, Mmn, Hqp);
  triplet_A.configure(opt);
  SingletOperator_BTDA_B singlet_B(epsilon_inv, Mmn, Hqp);
  singlet_B.configure(opt);
  BSE_OPERATOR<0, 0, 0, 1> triplet_B(epsilon_inv, Mmn, Hqp);
  triplet_B.configure(opt);

  ExcitonUKSOperator_TDA uks_A(epsilon_inv, Mspin, Hqp, Hqp);
  uks_A.configure(uopt);
  ExcitonUKSOperator_BTDA_B uks_B(epsilon_inv, Mspin, Hqp, Hqp);
  uks_B.configure(uopt);

  const Eigen::MatrixXd SA = Dense(singlet_A);
  const Eigen::MatrixXd TA = Dense(triplet_A);
  const Eigen::MatrixXd SB = Dense(singlet_B);
  const Eigen::MatrixXd TB = Dense(triplet_B);
  const Eigen::MatrixXd UA = Dense(uks_A);
  const Eigen::MatrixXd UB = Dense(uks_B);
  BOOST_REQUIRE_EQUAL(UA.rows(), 2 * SA.rows());
  BOOST_CHECK_SMALL((UA - UA.transpose()).cwiseAbs().maxCoeff(), 1e-12);
  BOOST_CHECK_SMALL((UB - UB.transpose()).cwiseAbs().maxCoeff(), 1e-12);
  BOOST_CHECK_SMALL((uks_A.diagonal() - UA.diagonal()).cwiseAbs().maxCoeff(),
                    1e-12);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_s(SA);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_t(TA);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_u(UA);
  const Eigen::VectorXd tda_ref =
      Sorted(Concat(es_s.eigenvalues(), es_t.eigenvalues()));
  const Eigen::VectorXd tda_uks = Sorted(es_u.eigenvalues());
  BOOST_CHECK_SMALL((tda_uks - tda_ref).cwiseAbs().maxCoeff(), 1e-10);

  const Eigen::VectorXd full_ref =
      Sorted(Concat(FullBSESquared(SA, SB), FullBSESquared(TA, TB)));
  const Eigen::VectorXd full_uks = Sorted(FullBSESquared(UA, UB));
  BOOST_CHECK_SMALL((full_uks - full_ref).cwiseAbs().maxCoeff(), 1e-10);
  libint2::finalize();
}

// The cached dense direct term must give the same products as the row loop.
BOOST_AUTO_TEST_CASE(direct_term_cache_equals_row_loop) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/bse/molecule.xyz");
  orbitals.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) + "/bse/3-21G.xml");
  AOBasis aobasis = orbitals.getDftBasis();
  orbitals.setNumberOfOccupiedLevels(4);
  Eigen::MatrixXd& MOs = orbitals.MOs().eigenvectors();
  MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/MOs.mm");
  const Eigen::MatrixXd Hqp = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse_operator/Hqp.mm");
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);
  Eigen::VectorXd epsilon_inv(aobasis.AOBasisSize());
  for (votca::Index i = 0; i < epsilon_inv.size(); ++i) {
    epsilon_inv(i) = 0.3 + 0.6 * double(i) / double(epsilon_inv.size());
  }
  BSEOperator_Options opt;
  opt.cmax = 8;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.rpamin = 0;
  opt.vmin = 0;

  std::srand(9);
  auto Compare = [&](auto& plain, auto& cached) {
    plain.configure(opt);
    cached.configure(opt);
    cached.set_direct_cache_limit(1e9);
    const Eigen::MatrixXd X = Eigen::MatrixXd::Random(plain.rows(), 5);
    const Eigen::MatrixXd ref = plain.matmul(X);
    const Eigen::MatrixXd first = cached.matmul(X);
    BOOST_CHECK(cached.direct_term_cached());
    BOOST_CHECK(!plain.direct_term_cached());
    const Eigen::MatrixXd second = cached.matmul(X);
    BOOST_CHECK_SMALL((first - ref).cwiseAbs().maxCoeff(), 1e-13);
    BOOST_CHECK_SMALL((second - ref).cwiseAbs().maxCoeff(), 1e-13);
  };
  SingletOperator_TDA s_plain(epsilon_inv, Mmn, Hqp);
  SingletOperator_TDA s_cached(epsilon_inv, Mmn, Hqp);
  Compare(s_plain, s_cached);
  SingletOperator_BTDA_B b_plain(epsilon_inv, Mmn, Hqp);
  SingletOperator_BTDA_B b_cached(epsilon_inv, Mmn, Hqp);
  Compare(b_plain, b_cached);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
