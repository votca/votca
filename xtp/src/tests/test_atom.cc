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

#define BOOST_TEST_MODULE atom_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/atom.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(atom_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  Atom atm(1, "C", Eigen::Vector3d::Zero());
}

BOOST_AUTO_TEST_CASE(getters_test) {
  Atom atm(3, "H", Eigen::Vector3d::Zero());
  BOOST_CHECK_EQUAL(atm.getId(), 3);
  BOOST_CHECK_EQUAL(atm.getName(), "H");
}

BOOST_AUTO_TEST_CASE(element_test) {
  Atom atm1(1, "CA3", Eigen::Vector3d::Zero());
  BOOST_CHECK_EQUAL(atm1.getElement(), "C");

  Atom atm2(2, "H", Eigen::Vector3d::Zero());
  BOOST_CHECK_EQUAL(atm2.getElement(), "H");

  BOOST_REQUIRE_THROW(Atom atm3(3, "Xacsdf", Eigen::Vector3d::Zero()),
                      std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
