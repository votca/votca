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

#define BOOST_TEST_MODULE eris_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/xtp/ERIs.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(eris_test)

BOOST_AUTO_TEST_CASE(fourcenter) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/eris/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/eris/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::MatrixXd dmat = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/dmat.mm");

  ERIs eris;
  eris.Initialize_4c(aobasis);
  Eigen::MatrixXd erissmall = eris.CalculateERIs_4c(dmat, 1e-20);

  Eigen::MatrixXd eris_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/eris_ref.mm");

  bool eris_check = erissmall.isApprox(eris_ref, 0.00001);
  if (!eris_check) {
    std::cout << "result eri" << std::endl;
    std::cout << erissmall << std::endl;
    std::cout << "ref eri" << std::endl;
    std::cout << eris_ref << std::endl;
    std::cout << " quotient" << std::endl;
    std::cout << erissmall.cwiseQuotient(eris_ref);
  }
  BOOST_CHECK_EQUAL(eris_check, 1);

  std::array<Eigen::MatrixXd, 2> both = eris.CalculateERIs_EXX_4c(dmat, 1e-20);
  const Eigen::MatrixXd& exx_small = both[1];
  Eigen::MatrixXd exx_ref = -votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/exx_ref.mm");

  Eigen::MatrixXd sum = both[0] + both[1];
  Eigen::MatrixXd sum_ref = exx_ref + eris_ref;
  bool sum_check = sum.isApprox(sum_ref, 1e-5);
  BOOST_CHECK_EQUAL(sum_check, 1);
  if (!sum_check) {
    std::cout << "result sum" << std::endl;
    std::cout << sum << std::endl;
    std::cout << "ref sum" << std::endl;
    std::cout << sum_ref << std::endl;
  }

  bool exxs_check = exx_small.isApprox(exx_ref, 0.00001);
  if (!eris_check) {
    std::cout << "result exx" << std::endl;
    std::cout << exx_small << std::endl;
    std::cout << "ref exx" << std::endl;
    std::cout << exx_ref << std::endl;
    std::cout << "quotient" << std::endl;
    std::cout << exx_small.cwiseQuotient(exx_ref);
  }
  BOOST_CHECK_EQUAL(exxs_check, 1);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(threecenter) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/eris/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/eris/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::MatrixXd dmat = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/dmat2.mm");

  Eigen::MatrixXd mos = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/mos.mm");

  ERIs eris;
  eris.Initialize(aobasis, aobasis);
  Eigen::MatrixXd exx_dmat =
      eris.CalculateERIs_EXX_3c(Eigen::MatrixXd::Zero(0, 0), dmat)[1];
  Eigen::MatrixXd exx_mo =
      eris.CalculateERIs_EXX_3c(mos.block(0, 0, 17, 4), dmat)[1];

  bool compare_exx = exx_mo.isApprox(exx_dmat, 1e-4);
  BOOST_CHECK_EQUAL(compare_exx, true);

  Eigen::MatrixXd exx_ref = -votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/exx_ref2.mm");

  bool compare_exx_ref = exx_ref.isApprox(exx_mo, 1e-5);
  if (!compare_exx_ref) {
    std::cout << "result exx" << std::endl;
    std::cout << exx_mo << std::endl;
    std::cout << "ref exx" << std::endl;
    std::cout << exx_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(compare_exx_ref, true);

  Eigen::MatrixXd eri = eris.CalculateERIs_3c(dmat);

  Eigen::MatrixXd eris_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/eris_ref2.mm");

  bool compare_eris = eris_ref.isApprox(eri, 1e-5);
  if (!compare_eris) {
    std::cout << "result eris" << std::endl;
    std::cout << eri << std::endl;
    std::cout << "ref eris" << std::endl;
    std::cout << eris_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(compare_eris, true);
}

BOOST_AUTO_TEST_SUITE_END()
