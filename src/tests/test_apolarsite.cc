/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE apolararsite_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/apolarsite.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(apolararsite_test)

BOOST_AUTO_TEST_CASE(constructors_test) { APolarSite ps(1, "ps1"); }

BOOST_AUTO_TEST_CASE(getters_test) {
  APolarSite ps(1, "ps2");
  BOOST_CHECK_EQUAL(ps.getId(), 1);
  BOOST_CHECK_EQUAL(ps.getName(), "ps2");
}

BOOST_AUTO_TEST_SUITE_END()
