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
#define BOOST_TEST_MODULE settings_test

// Standard includes
#include <fstream>
#include <iostream>
#include <string>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/settings.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(settings_test)

BOOST_AUTO_TEST_CASE(create_settings) {

  Settings qmpackage_template{"package"};
  qmpackage_template.load_from_xml(std::string(XTP_TEST_DATA_FOLDER) +
                                   "/settings/defaults.xml");
  auto basisset = qmpackage_template.get("basisset");
  auto orca_guess = qmpackage_template.get("orca.scf");
  BOOST_TEST(basisset == "ubecppol");
  BOOST_TEST(orca_guess == "GUESS PMODEL");
}

BOOST_AUTO_TEST_CASE(test_amend) {

  Settings qmpackage_template{"package"};
  qmpackage_template.load_from_xml(std::string(XTP_TEST_DATA_FOLDER) +
                                   "/settings/defaults2.xml");
  Settings user_input("package");
  user_input.load_from_xml(std::string(XTP_TEST_DATA_FOLDER) +
                           "/settings/user_input.xml");
  user_input.amend(qmpackage_template);
  user_input.add("orca.property", "42");
  user_input.validate();

  auto basisset = user_input.get("functional");
  auto ecp = user_input.get("ecp");
  auto orca_guess = user_input.get("orca.scf");
  auto executable = user_input.get("executable");
  auto orca_prop = user_input.get("orca.property");
  BOOST_TEST(basisset == "B3LYP");
  BOOST_TEST(ecp == "ecp");
  BOOST_TEST(executable == "/opt/orca/orca");
  BOOST_TEST(orca_guess == "GUESS PMODEL");
  BOOST_TEST(orca_prop == "42");
}

BOOST_AUTO_TEST_SUITE_END()
