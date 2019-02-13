/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#define BOOST_TEST_MODULE vec_test
#include <boost/test/unit_test.hpp>
#include <votca/tools/vec.h>

using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(vec_test)

BOOST_AUTO_TEST_CASE(overloadoperator_test) {

  vec v1(1, 1, 1);
  vec v2(1.0);
  vec v3(1.1);

  BOOST_CHECK(v1.isClose(v2, 0.001));

  BOOST_CHECK(!v1.isClose(v3, 0.001));
}

BOOST_AUTO_TEST_CASE(Eigenconv_test) {

  vec v1(1, 0, 0);
  Eigen::Vector3d unit = Eigen::Vector3d::UnitX();

  Eigen::Vector3d conv = vec(unit).toEigen();

  BOOST_CHECK(v1.toEigen().isApprox(unit, 0.001));

  BOOST_CHECK(conv.isApprox(unit, 0.0001));
}

BOOST_AUTO_TEST_SUITE_END()
