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

#define BOOST_TEST_MODULE sigma_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigma_exact.h"
#include "votca/xtp/threecenter.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(sigma_test)

BOOST_AUTO_TEST_CASE(sigma_full) {

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << 0.0468207, 0.0907801, 0.0907801, 0.104563, 0.592491, 0.663355,
      0.663355, 0.768373, 1.69292, 1.97724, 1.97724, 2.50877, 2.98732, 3.4418,
      3.4418, 4.81084, 17.1838;

  Eigen::MatrixXd MOs =  votca::tools::EigenIO_MatrixMarket::ReadMatrix(
  std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/MOs.mm");

  Logger log;
  TCMatrix_gwbse Mmn{log};
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);

  RPA rpa(log, Mmn);
  rpa.setRPAInputEnergies(mo_energy);
  rpa.configure(4, 0, 16);

  Sigma_Exact sigma = Sigma_Exact(Mmn, rpa);
  Sigma_Exact::options opt;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.qpmax = 16;
  opt.rpamin = 0;
  opt.rpamax = 16;
  sigma.configure(opt);

  Eigen::MatrixXd x = sigma.CalcExchangeMatrix();

  Eigen::MatrixXd x_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
  std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/x_ref.mm");

  bool check_x = x_ref.isApprox(x, 1e-5);
  if (!check_x) {
    cout << "Sigma X" << endl;
    cout << x << endl;
    cout << "Sigma X ref" << endl;
    cout << x_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_x, true);

  sigma.PrepareScreening();
  Eigen::MatrixXd c = sigma.CalcCorrelationOffDiag(mo_energy);
  c.diagonal() = sigma.CalcCorrelationDiag(mo_energy);

  Eigen::MatrixXd c_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
  std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/c_ref.mm");

  bool check_c_diag = c.diagonal().isApprox(c_ref.diagonal(), 1e-5);
  if (!check_c_diag) {
    cout << "Sigma C" << endl;
    cout << c.diagonal() << endl;
    cout << "Sigma C ref" << endl;
    cout << c_ref.diagonal() << endl;
  }
  BOOST_CHECK_EQUAL(check_c_diag, true);

  bool check_c = c.isApprox(c_ref, 1e-5);
  if (!check_c) {
    cout << "Sigma C" << endl;
    cout << c << endl;
    cout << "Sigma C ref" << endl;
    cout << c_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_c, true);
}

BOOST_AUTO_TEST_SUITE_END()
