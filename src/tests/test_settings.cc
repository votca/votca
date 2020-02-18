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
#define BOOST_TEST_MODULE settings_test

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <string>
#include <votca/xtp/settings.h>
using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(settings_test)

BOOST_AUTO_TEST_CASE(create_settings) {
  Settings qmpackage_template{"package"};
  std::string env = getenv("VOTCASHARE");
  auto xmlFile = std::string(getenv("VOTCASHARE")) +
                 std::string("/xtp/packages/qmpackage_template.xml");
  qmpackage_template.load_from_xml(xmlFile);

  auto gaussian_memory = qmpackage_template.get("gaussian.memory");
  auto basisset = qmpackage_template.get("basisset");
  BOOST_TEST(gaussian_memory == "1GB");
  BOOST_TEST(basisset == "ubecppol");
}
BOOST_AUTO_TEST_SUITE_END()
