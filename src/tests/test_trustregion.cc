/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE trustregion_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/xtp/trustregion.h>

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(trustregion_test)

BOOST_AUTO_TEST_CASE(subproblem_test) {

  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(5);
  gradient << 5.92833, 9.88054, 9.88054, 9.88054, 9.88054;
  Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(5, 5);
  hessian << 1.0835, 0.139166, 0.139166, 0.139166, 0.139166, 0.139166, 1.23194,
      0.231943, 0.231943, 0.231943, 0.139166, 0.231943, 1.23194, 0.231943,
      0.231943, 0.139166, 0.231943, 0.231943, 1.23194, 0.231943, 0.139166,
      0.231943, 0.231943, 0.231943, 1.23194;

  double radius = 0.1;
  TrustRegion trf;
  Eigen::VectorXd newstep = trf.CalculateStep(gradient, hessian, radius);

  bool equal = std::abs(newstep.norm() - radius) < 1e-9;
  if (!equal) {
    std::cout << "newstep_norm:" << newstep.norm() << " trust_radius:" << radius
              << std::endl;
  }

  BOOST_CHECK_EQUAL(equal, 1);
}

BOOST_AUTO_TEST_SUITE_END()
