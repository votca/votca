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

#define BOOST_TEST_MODULE adiis_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/adiis.h"
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(adiis_test)

BOOST_AUTO_TEST_CASE(coeffs_test) {
  std::vector<Eigen::MatrixXd> _dmatHist;
  std::vector<Eigen::MatrixXd> _mathist;

  ADIIS adiis;

  for (votca::Index i = 0; i < 3; i++) {
    _dmatHist.push_back(votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/adiis/dmatHist_" +
        std::to_string(i) + ".mm"));
    _mathist.push_back(votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/adiis/mathist_" +
        std::to_string(i) + ".mm"));
  }

  Eigen::VectorXd Coeffs = adiis.CalcCoeff(_dmatHist, _mathist);

  Eigen::VectorXd Ref = Eigen::VectorXd::Zero(3);
  Ref << 4.45639501e-17, 1.76102089e-18, 1;

  bool check_adiis = Coeffs.isApprox(Ref, 0.00001);
  if (!check_adiis) {
    std::cout << "Ref:" << Ref << std::endl;
    std::cout << "Coeffs:" << Coeffs << std::endl;
  }

  BOOST_CHECK_EQUAL(check_adiis, 1);
}

BOOST_AUTO_TEST_SUITE_END()
