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

#define BOOST_TEST_MODULE rpa_test

// Third party includes
#include "boost/test/unit_test.hpp"

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/threecenter.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(rpa_test)

BOOST_AUTO_TEST_CASE(rpa_calcenergies) {

  Logger log;
  TCMatrix_gwbse Mmn;
  Eigen::VectorXd eigenvals;
  RPA rpa(log, Mmn);
  rpa.configure(4, 0, 9);
  Eigen::VectorXd dftenergies = Eigen::VectorXd::Zero(10);
  dftenergies << -0.5, -0.4, -0.3, -0.2, -0.2, -0.1, 0, 0.1, 0.2, 0.3;
  Eigen::VectorXd gwenergies = Eigen::VectorXd::Zero(7);
  gwenergies << -0.15, -0.05, 0.05, 0.15, 0.45, 0.55, 0.65;
  votca::Index qpmin = 1;
  rpa.UpdateRPAInputEnergies(dftenergies, gwenergies, qpmin);
  Eigen::VectorXd rpaenergies = rpa.getRPAInputEnergies();
  Eigen::VectorXd rpaenergies_ref = Eigen::VectorXd::Zero(10);
  rpaenergies_ref << -0.85, -0.15, -0.05, 0.05, 0.15, 0.45, 0.55, 0.65, 0.75,
      0.85;
  bool e_check = rpaenergies_ref.isApprox(rpaenergies, 0.0001);

  if (!e_check) {
    cout << "energy" << endl;
    cout << rpaenergies << endl;
    cout << "energy_ref" << endl;
    cout << rpaenergies_ref << endl;
  }
  BOOST_CHECK_EQUAL(e_check, true);
}

BOOST_AUTO_TEST_CASE(rpa_full) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/rpa/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/rpa/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::VectorXd eigenvals = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/eigenvals.mm");

  Eigen::MatrixXd eigenvectors = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/eigenvectors.mm");
  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, eigenvectors);

  RPA rpa(log, Mmn);
  rpa.configure(4, 0, 16);
  rpa.setRPAInputEnergies(eigenvals);
  Eigen::MatrixXd e_i = rpa.calculate_epsilon_i(0.5);

  Eigen::MatrixXd i_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/i_ref.mm");
  bool i_check = i_ref.isApprox(e_i, 0.0001);

  if (!i_check) {
    cout << "Epsilon_i" << endl;
    cout << e_i << endl;
    cout << "Epsilon_i_ref" << endl;
    cout << i_ref << endl;
  }
  BOOST_CHECK_EQUAL(i_check, 1);

  Eigen::MatrixXd e_r = rpa.calculate_epsilon_r(0.0);

  Eigen::MatrixXd r_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/r_ref.mm");
  bool r_check = r_ref.isApprox(e_r, 0.0001);

  if (!r_check) {
    cout << "Epsilon_r" << endl;
    cout << e_r << endl;
    cout << "Epsilon_r_ref" << endl;
    cout << r_ref << endl;
  }

  BOOST_CHECK_EQUAL(r_check, 1);

  Eigen::MatrixXd e_r_complex =
      rpa.calculate_epsilon_r(std::complex<double>(0.5, 0.5));

  Eigen::MatrixXd r_complex_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/rpa/r_complex_ref.mm");
  bool r_complex_check = r_complex_ref.isApprox(e_r_complex, 0.0001);

  if (!r_complex_check) {
    cout << "Epsilon_r_complex" << endl;
    cout << e_r_complex << endl;
    cout << "Epsilon_r_compelx_ref" << endl;
    cout << r_complex_ref << endl;
  }

  BOOST_CHECK_EQUAL(r_complex_check, 1);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
