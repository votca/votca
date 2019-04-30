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
#define BOOST_TEST_MODULE parameters_test

#include "../../include/votca/tools/parameters.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(parameters_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Parameters parameters; }

BOOST_AUTO_TEST_CASE(set_test) {
  Parameters parameters;

  double mass = 2.0;
  parameters.set(Parameters::Parameter::Mass, mass);
}

BOOST_AUTO_TEST_CASE(get_test) {
  Parameters parameters;

  double mass = 2.0;
  parameters.set(Parameters::Parameter::Mass, mass);
  double mass_check = parameters.get<double>(Parameters::Parameter::Mass);
  BOOST_CHECK_EQUAL(mass, mass_check);
}
/*
BOOST_AUTO_TEST_CASE(export_test) {
  int num = 87094;
  string name = "John Doe";
  double height = 167.8;
  string address = "12 Koogler St, Dayton, OH 32345";
  int age = 24;
  string fav_col = "mauve";
  Person John(num, name, address, age, fav_col, height);

  TypeConverter<Person, ContactInfo> converter;

  converter.importData(John);
  ContactInfo cont_info;
  converter.exportData(cont_info);

  BOOST_CHECK_EQUAL(cont_info.getPhoneNum(), num);
  BOOST_CHECK_EQUAL(cont_info.getName(), name);
  BOOST_CHECK_EQUAL(cont_info.getAddress(), address);
}*/

BOOST_AUTO_TEST_SUITE_END()
