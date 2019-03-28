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

#define BOOST_TEST_MODULE property_test
#include <boost/test/unit_test.hpp>
#include <votca/tools/property.h>

using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(property_test)

BOOST_AUTO_TEST_CASE(eigen_test) {
  Property prop;
  prop.add("vec", "1 2 3 4 5 6 7 8");
  Eigen::VectorXd result = prop.get("vec").as<Eigen::VectorXd>();

  Eigen::VectorXd ref;
  ref.resize(8);
  ref << 1, 2, 3, 4, 5, 6, 7, 8;
  BOOST_CHECK_EQUAL(ref.isApprox(result, 0.0001), true);

  BOOST_CHECK_THROW(prop.get("vec").as<Eigen::Vector3d>();, runtime_error);

  Property prop2;
  prop.add("vec", "1 2 3");
  Eigen::Vector3d result2 = prop.get("vec").as<Eigen::Vector3d>();
  Eigen::Vector3d ref2{1, 2, 3};
  BOOST_CHECK_EQUAL(ref2.isApprox(result2, 0.0001), true);
}

BOOST_AUTO_TEST_SUITE_END()
