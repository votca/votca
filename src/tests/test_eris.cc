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

BOOST_AUTO_TEST_CASE(fourcenter_cache) {

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
  eris.Initialize_4c_small_molecule(aobasis);
  Mat_p_Energy erissmall = eris.CalculateERIs_4c_small_molecule(dmat);

  Eigen::MatrixXd eris_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/eris_ref.mm");

  bool eris_check = erissmall.matrix().isApprox(eris_ref, 0.00001);
  BOOST_CHECK_EQUAL(eris_check, 1);

  Mat_p_Energy exx_small = eris.CalculateEXX_4c_small_molecule(dmat);

  Eigen::MatrixXd exx_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/exx_ref.mm");

  bool exxs_check = exx_small.matrix().isApprox(exx_ref, 0.00001);
  BOOST_CHECK_EQUAL(exxs_check, 1);
}

BOOST_AUTO_TEST_CASE(threecenter) {

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
  Mat_p_Energy exx_dmat = eris.CalculateEXX(dmat);
  Mat_p_Energy exx_mo = eris.CalculateEXX(mos.block(0, 0, 17, 4), dmat);

  bool compare_exx = exx_mo.matrix().isApprox(exx_dmat.matrix(), 1e-4);
  BOOST_CHECK_EQUAL(compare_exx, true);

  Eigen::MatrixXd exx_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/exx_ref2.mm");

  bool compare_exx_ref = exx_ref.isApprox(exx_mo.matrix(), 1e-5);
  if (!compare_exx_ref) {
    std::cout << "result exx" << std::endl;
    std::cout << exx_mo.matrix() << std::endl;
    std::cout << "ref exx" << std::endl;
    std::cout << exx_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(compare_exx_ref, true);

  Mat_p_Energy eri = eris.CalculateERIs(dmat);

  Eigen::MatrixXd eris_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/eris_ref2.mm");

  bool compare_eris = eris_ref.isApprox(eri.matrix(), 1e-5);
  if (!compare_eris) {
    std::cout << "result eris" << std::endl;
    std::cout << eri.matrix() << std::endl;
    std::cout << "ref eris" << std::endl;
    std::cout << eris_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(compare_eris, true);
}

BOOST_AUTO_TEST_CASE(fourcenter_direct) {

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/eris/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/eris/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::MatrixXd dmat = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/eris/dmat3.mm");

  ERIs eris1;
  ERIs eris2;
  eris1.Initialize_4c_screening(aobasis, 1e-10);
  eris2.Initialize_4c_small_molecule(aobasis);
  Mat_p_Energy eris_direct = eris1.CalculateERIs_4c_direct(aobasis, dmat);
  Mat_p_Energy eris_cached = eris2.CalculateERIs_4c_small_molecule(dmat);

  bool check_eris = eris_direct.matrix().isApprox(eris_cached.matrix(), 0.001);
  if (!check_eris) {
    std::cout << "eris_direct.matrix()" << std::endl;
    std::cout << eris_direct.matrix() << std::endl;
    std::cout << "eris_cached.matrix()" << std::endl;
    std::cout << eris_cached.matrix() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_eris, 1);
}

BOOST_AUTO_TEST_SUITE_END()
