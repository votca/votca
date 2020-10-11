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

#define BOOST_TEST_MODULE threecenter_gwbse_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA inlcudes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/threecenter.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(threecenter_gwbse_test)
BOOST_AUTO_TEST_CASE(threecenter_gwbse) {
  libint2::initialize();
  QMMolecule mol(" ", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/threecenter_gwbse/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) +
             "/threecenter_gwbse/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);

  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/MOs.mm");

  Logger log;
  TCMatrix_gwbse tc;
  tc.Initialize(aobasis.AOBasisSize(), 0, 5, 0, 7);
  tc.Fill(aobasis, aobasis, MOs);

  Eigen::MatrixXd ref0b = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/ref0b.mm");

  bool check0_before = ref0b.isApprox(tc[0], 1e-5);
  if (!check0_before) {
    cout << "tc0" << endl;
    cout << tc[0] << endl;
    cout << "tc0_ref" << endl;
    cout << ref0b << endl;
  }
  BOOST_CHECK_EQUAL(check0_before, true);

  Eigen::MatrixXd ref2b = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/ref2b.mm");

  bool check2_before = ref2b.isApprox(tc[2], 1e-5);
  if (!check2_before) {
    cout << "tc2" << endl;
    cout << tc[2] << endl;
    cout << "tc2_ref" << endl;
    cout << ref2b << endl;
  }

  BOOST_CHECK_EQUAL(check2_before, true);

  Eigen::MatrixXd ref4b = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_gwbse/ref4b.mm");

  bool check4_before = ref4b.isApprox(tc[4], 1e-5);
  if (!check4_before) {
    cout << "tc4" << endl;
    cout << tc[4] << endl;
    cout << "tc4_ref" << endl;
    cout << ref4b << endl;
  }

  BOOST_CHECK_EQUAL(check4_before, true);

  Eigen::MatrixXd auxmatrix =
      Eigen::MatrixXd::Identity(aobasis.AOBasisSize(), aobasis.AOBasisSize());
  tc.MultiplyRightWithAuxMatrix(auxmatrix);

  bool check0_after = ref0b.isApprox(tc[0], 1e-5);
  if (!check0_after) {
    cout << "tc0" << endl;
    cout << tc[0] << endl;
    cout << "tc0_ref" << endl;
    cout << ref0b << endl;
  }
  BOOST_CHECK_EQUAL(check0_after, true);

  bool check2_after = ref2b.isApprox(tc[2], 1e-5);
  if (!check2_after) {
    cout << "tc2" << endl;
    cout << tc[2] << endl;
    cout << "tc2_ref" << endl;
    cout << ref2b << endl;
  }

  BOOST_CHECK_EQUAL(check2_after, true);

  bool check4_after = ref4b.isApprox(tc[4], 1e-5);
  if (!check4_after) {
    cout << "tc4" << endl;
    cout << tc[4] << endl;
    cout << "tc4_ref" << endl;
    cout << ref4b << endl;
  }

  BOOST_CHECK_EQUAL(check4_after, true);

  libint2::finalize();
}
BOOST_AUTO_TEST_SUITE_END()
