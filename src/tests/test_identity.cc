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

#define BOOST_TEST_MODULE identity_test
#include <boost/test/unit_test.hpp>
#include <exception>
#include <votca/tools/identity.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(identity_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  Identity<int> id;
  Identity<long int> id2(232);
}

BOOST_AUTO_TEST_CASE(simple_test) {
  Identity<int> id;
  Identity<int> id2(32);
  BOOST_CHECK_EQUAL(id2.getId(), 32);
  id2.setId(34);
  BOOST_CHECK_EQUAL(id2.getId(), 34);
}

BOOST_AUTO_TEST_SUITE_END()
