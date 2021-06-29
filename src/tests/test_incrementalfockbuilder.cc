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
#include <boost/test/tools/old/interface.hpp>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE incrementalfockbuilder_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/IncrementalFockBuilder.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(incrementalfockbuilder)

BOOST_AUTO_TEST_CASE(off_test) {
  Logger log;

  // switched off
  IncrementalFockBuilder fb(log, 0.0, 3);

  Eigen::MatrixXd dmat = Eigen::MatrixXd::Random(10, 10);
  fb.Configure(dmat);
  fb.Start(2, 1e-5);

  Eigen::MatrixXd J = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd K = Eigen::MatrixXd::Random(10, 10);
  fb.resetMatrices(J, K, dmat);
  BOOST_CHECK(J.isApproxToConstant(0.0));
  BOOST_CHECK(K.isApproxToConstant(0.0));

  BOOST_CHECK(fb.getDmat_diff().isApprox(dmat));

  fb.UpdateDmats(dmat, 1e-5, 2);

  // basically checks that using the same dmat twice leads to a diff of zero
  BOOST_CHECK(fb.getDmat_diff().isApproxToConstant(0.0));
}

BOOST_AUTO_TEST_CASE(on_test) {
  Logger log;

  // switched on
  IncrementalFockBuilder fb(log, 1e-5, 2);

  Eigen::MatrixXd dmat = Eigen::MatrixXd::Random(10, 10);
  fb.Configure(dmat);
  votca::Index iteration = 2;
  fb.Start(iteration, 1e-6);

  Eigen::MatrixXd J = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd J2 = J;
  Eigen::MatrixXd K = Eigen::MatrixXd::Random(10, 10);
  Eigen::MatrixXd K2 = K;

  fb.resetMatrices(J, K, dmat);
  BOOST_CHECK(J.isApprox(J2));
  BOOST_CHECK(K.isApprox(K2));

  fb.UpdateCriteria(1e-6, iteration);
  fb.UpdateDmats(dmat, 1e-6, iteration);
  // basically checks that using the same dmat twice leads to a diff of zero
  BOOST_CHECK(fb.getDmat_diff().isApproxToConstant(0.0));
  iteration++;
  fb.UpdateCriteria(1e-6, iteration);
  fb.UpdateDmats(dmat, 1e-6, iteration);
  iteration++;
  fb.UpdateCriteria(1e-6, iteration);
  fb.UpdateDmats(dmat, 1e-6, iteration);

  // now reset should trigger
  fb.resetMatrices(J, K, dmat);
  BOOST_CHECK(J.isApproxToConstant(0.0));
  BOOST_CHECK(K.isApproxToConstant(0.0));
}

BOOST_AUTO_TEST_SUITE_END()
