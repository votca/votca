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

#define BOOST_TEST_MODULE anderson_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/anderson_mixing.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(anderson_test)

BOOST_AUTO_TEST_CASE(coeffs_test) {

  Anderson _mixing;
  _mixing.Configure(3, 0.7);
  Eigen::VectorXd in1 = Eigen::VectorXd::Zero(7);
  in1 << -0.580533, -0.535803, -0.476481, -0.380558, 0.0969526, 0.133036,
      0.164243;
  _mixing.UpdateInput(in1);

  Eigen::VectorXd out1 = Eigen::VectorXd::Zero(7);
  out1 << -0.604342, -0.548675, -0.488088, -0.385654, 0.106193, 0.139172,
      0.170433;
  _mixing.UpdateOutput(out1);

  Eigen::VectorXd mixed = _mixing.MixHistory();

  Eigen::VectorXd ref1 = Eigen::VectorXd::Zero(7);
  ref1 << -0.597199, -0.544813, -0.484606, -0.384126, 0.103421, 0.137331,
      0.168576;

  bool check_linear = mixed.isApprox(ref1, 0.00001);
  if (!check_linear) {
    std::cout << "Ref:" << ref1 << std::endl;
    std::cout << "Linear:" << mixed << std::endl;
  }

  BOOST_CHECK_EQUAL(check_linear, 1);

  _mixing.UpdateInput(mixed);

  Eigen::VectorXd out2 = Eigen::VectorXd::Zero(7);
  out2 << -0.605576, -0.549458, -0.488876, -0.385821, 0.106788, 0.139509,
      0.170768;
  _mixing.UpdateOutput(out2);

  mixed = _mixing.MixHistory();

  Eigen::VectorXd ref2 = Eigen::VectorXd::Zero(7);
  ref2 << -0.606303, -0.549862, -0.489247, -0.385968, 0.10708, 0.139698,
      0.170959;

  bool check_nonlinear_order2 = mixed.isApprox(ref2, 0.00001);
  if (!check_nonlinear_order2) {
    std::cout << "Ref:" << ref2 << std::endl;
    std::cout << "Nonlinear 2nd order:" << mixed << std::endl;
  }

  BOOST_CHECK_EQUAL(check_nonlinear_order2, 1);
  _mixing.UpdateInput(mixed);

  Eigen::VectorXd out3 = Eigen::VectorXd::Zero(7);
  out3 << -0.606242, -0.549887, -0.489296, -0.385898, 0.107162, 0.139718,
      0.170977;
  _mixing.UpdateOutput(out3);

  mixed = _mixing.MixHistory();
  Eigen::VectorXd ref3 = Eigen::VectorXd::Zero(7);
  ref3 << -0.606242, -0.549888, -0.489296, -0.385897, 0.107163, 0.139718,
      0.170977;
  bool check_nonlinear_order3 = mixed.isApprox(ref3, 0.00001);
  if (!check_nonlinear_order3) {
    std::cout << "Ref:" << ref3 << std::endl;
    std::cout << "Nonlinear 3rd order:" << mixed << std::endl;
  }

  BOOST_CHECK_EQUAL(check_nonlinear_order3, 1);
}

BOOST_AUTO_TEST_SUITE_END()
