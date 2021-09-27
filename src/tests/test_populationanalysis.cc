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
#include <libint2/initialize.h>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE populationanalysis_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/populationanalysis.h"

using namespace std;
using namespace votca::xtp;
using namespace votca;
BOOST_AUTO_TEST_SUITE(populationanalysis_test)

BOOST_AUTO_TEST_CASE(atompop) {
  libint2::initialize();
  Orbitals orb;
  orb.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                             "/populationanalysis/molecule.xyz");
  orb.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                      "/populationanalysis/3-21G.xml");
  orb.setBasisSetSize(17);
  orb.setNumberOfOccupiedLevels(5);

  Eigen::MatrixXd& MOs = orb.MOs().eigenvectors();
  orb.MOs().eigenvalues() = Eigen::VectorXd::Ones(17);
  MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/populationanalysis/MOs.mm");

  Orbitals orb2 = orb;
  QMState s("n");
  Lowdin low;
  StaticSegment result = low.CalcChargeperAtom(orb, s);

  Eigen::VectorXd charge = Eigen::VectorXd::Zero(result.size());
  for (votca::Index i = 0; i < votca::Index(result.size()); i++) {
    charge(i) = result[i].getCharge();
  }

  Eigen::VectorXd charge_ref = Eigen::VectorXd::Zero(5);
  charge_ref << 0.68862, -0.172155, -0.172154, -0.172154, -0.172155;

  bool check_lowdin = charge_ref.isApprox(charge, 1e-5);
  BOOST_CHECK_EQUAL(check_lowdin, true);
  if (!check_lowdin) {
    cout << "charge" << endl;
    cout << charge << endl;
    cout << "chargeref" << endl;
    cout << charge_ref << endl;
  }

  Mulliken mul;
  StaticSegment result2 = mul.CalcChargeperAtom(orb2, s);
  Eigen::VectorXd charge2 = Eigen::VectorXd::Zero(result2.size());
  for (votca::Index i = 0; i < votca::Index(result2.size()); i++) {
    charge2(i) = result2[i].getCharge();
  }

  Eigen::VectorXd charge_ref2 = Eigen::VectorXd::Zero(5);
  charge_ref2 << 1.21228, -0.303069, -0.303067, -0.303068, -0.303068;

  bool check_mulliken = charge_ref2.isApprox(charge2, 1e-5);
  BOOST_CHECK_EQUAL(check_mulliken, true);
  if (!check_mulliken) {
    cout << "charge" << endl;
    cout << charge2 << endl;
    cout << "chargeref" << endl;
    cout << charge_ref2 << endl;
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(fragment_pop) {
  libint2::initialize();
  Orbitals orb;
  orb.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                             "/populationanalysis/molecule.xyz");
  orb.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                      "/populationanalysis/3-21G.xml");
  orb.setBasisSetSize(17);
  orb.setNumberOfOccupiedLevels(5);
  orb.MOs().eigenvalues() = Eigen::VectorXd::Ones(17);

  Eigen::MatrixXd MOs2 = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/populationanalysis/MOs2.mm");

  orb.MOs().eigenvectors() = MOs2;
  Eigen::VectorXd se_ref = Eigen::VectorXd::Zero(3);
  se_ref << 0.107455, 0.107455, 0.107455;
  orb.BSESinglets().eigenvalues() = se_ref;

  // reference coefficients
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/populationanalysis/spsi_ref.mm");

  orb.BSESinglets().eigenvectors() = spsi_ref;
  orb.setBSEindices(0, 16);
  orb.setTDAApprox(true);

  std::vector<QMFragment<BSE_Population> > frags;
  QMFragment<BSE_Population> f1(0, "0 1");
  QMFragment<BSE_Population> f2(1, "2 3 4");
  frags.push_back(f1);
  frags.push_back(f2);

  orb.DensityMatrixGroundState();
  QMState n("n");
  Lowdin low;
  BOOST_REQUIRE_THROW(low.CalcChargeperFragment(frags, orb, n.Type()),
                      std::runtime_error);
  QMState s1("s1");
  low.CalcChargeperFragment(frags, orb, s1.Type());

  BOOST_CHECK_CLOSE(frags[0].value().Gs, 0.5164649, 1e-5);
  BOOST_CHECK_CLOSE(frags[1].value().Gs, -0.5164628, 1e-5);

  Eigen::VectorXd f1E_ref = Eigen::VectorXd::Zero(3);
  f1E_ref << -0.384176, -0.812396, -0.414518;

  Eigen::VectorXd f1H_ref = Eigen::VectorXd::Zero(3);
  f1H_ref << 0.657215, 0.622434, 0.654751;

  Eigen::VectorXd f2E_ref = Eigen::VectorXd::Zero(3);
  f2E_ref << -0.615827, -0.187602, -0.58548;

  Eigen::VectorXd f2H_ref = Eigen::VectorXd::Zero(3);
  f2H_ref << 0.342785, 0.377565, 0.34525;

  bool check_f1e = frags[0].value().E.isApprox(f1E_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_f1e, true);
  if (!check_f1e) {
    cout << "charge" << endl;
    cout << frags[0].value().E << endl;
    cout << "chargeref" << endl;
    cout << f1E_ref << endl;
  }

  bool check_f1h = frags[0].value().H.isApprox(f1H_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_f1h, true);
  if (!check_f1h) {
    cout << "charge" << endl;
    cout << frags[0].value().H << endl;
    cout << "chargeref" << endl;
    cout << f1H_ref << endl;
  }
  bool check_f2e = frags[1].value().E.isApprox(f2E_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_f2e, true);
  if (!check_f2e) {
    cout << "charge" << endl;
    cout << frags[1].value().E << endl;
    cout << "chargeref" << endl;
    cout << f2E_ref << endl;
  }
  bool check_f2h = frags[1].value().H.isApprox(f2H_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_f2h, true);
  if (!check_f2h) {
    cout << "charge" << endl;
    cout << frags[1].value().H << endl;
    cout << "chargeref" << endl;
    cout << f2H_ref << endl;
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
