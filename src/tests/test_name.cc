/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE name_test
#include <boost/test/unit_test.hpp>
#include <exception>
#include <string>
#include <votca/tools/name.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(name_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  Name nm;
  Name nm2("S");
}

BOOST_AUTO_TEST_CASE(accessors_test) {
  Name nm;
  BOOST_CHECK_THROW(nm.getName(), runtime_error);
  nm.setName("New Name");
  BOOST_CHECK_EQUAL(nm.getName(), "New Name");
  Name nm2("Name2");
  BOOST_CHECK_EQUAL(nm2.getName(), "Name2");
}

BOOST_AUTO_TEST_SUITE_END()
