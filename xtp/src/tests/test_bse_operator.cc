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
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/bse_operator.h"
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

BOOST_AUTO_TEST_SUITE_END()
