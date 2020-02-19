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
#include <fstream>
#include <iostream>
#include <string>
#include <votca/xtp/settings.h>
using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(settings_test)

BOOST_AUTO_TEST_CASE(create_settings) {

  std::ofstream defaults("defaults.xml");

  defaults << "<package>\n"
           << "<name>orca</name>\n"
           << "<charge>0</charge>\n"
           << "<spin>1</spin>\n"
           << "<basisset>ubecppol</basisset>\n"
           << "<auxbasisset>aux-ubecppol</auxbasisset>\n"
           << "<optimize>false</optimize>\n"
           << "<functional>pbe</functional>\n"
           << "<scratch>/tmp/xtp_qmpackage</scratch>\n"
           << "<polarisation>false</polarisation>\n"
           << "<gaussian>\n"
           << "    <memory>1GB</memory>\n"
           << "</gaussian>\n"
           << "<orca>\n"
           << "   <scf>GUESS PMODEL</scf>\n"
           << "</orca>\n"
           << "</package>";
  defaults.close();

  Settings qmpackage_template{"package"};
  qmpackage_template.load_from_xml("defaults.xml");
  auto gaussian_memory = qmpackage_template.get("gaussian.memory");
  auto basisset = qmpackage_template.get("basisset");
  auto orca_guess = qmpackage_template.get("orca.scf");
  BOOST_TEST(gaussian_memory == "1GB");
  BOOST_TEST(basisset == "ubecppol");
  BOOST_TEST(orca_guess == "GUESS PMODEL");
}

BOOST_AUTO_TEST_CASE(create_section) {
  std::ofstream orca_prop("orca_prop.xml");
  orca_prop << "<package>\n"
            << "<orca>\n"
            << "   <scf>GUESS PMODEL</scf>\n"
            << "</orca>\n"
            << "</package>\n";
  orca_prop.close();

  Settings orca_settings{"package"};
  orca_settings.load_from_xml("orca_prop.xml");
  std::string section = orca_settings.CreateInputSection("orca.scf");
  std::cout << "section:\n" << section << "\n";
}

BOOST_AUTO_TEST_SUITE_END()
