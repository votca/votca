/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE calculator_test
#include "../../include/votca/tools/calculator.h"
#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>
#include <fstream>

using namespace ::votca;

BOOST_AUTO_TEST_SUITE(calculator_test)

BOOST_AUTO_TEST_CASE(load_defaults_test) {

  class TestCalc : public tools::Calculator {

   public:
    std::string Identify() override { return "testcalc"; }

    void Initialize(const tools::Property &user_options) override {

      // Create folder for test
      const char dir_path[] = "calculators";
      boost::filesystem::path dir(dir_path);
      boost::filesystem::create_directory(dir);
      dir.append("xml");
      boost::filesystem::create_directory(dir);

      std::ofstream defaults("calculators/xml/testcalc.xml");
      defaults << "<options>\n"
               << "<testcalc>\n"
               << "<option1>0</option1>\n"
               << "<option2>3.141592</option2>\n"
               << "</testcalc>\n"
               << "</options>";
      defaults.close();

      // Load and check the options
      tools::Property final_opt =
          LoadDefaultsAndUpdateWithUserOptions("calculators", user_options);

      Index prop1 = final_opt.get("option1").as<votca::Index>();
      std::string prop2 = final_opt.get("option2").as<std::string>();
      std::string prop3 = final_opt.get("option3.nested").as<std::string>();
      BOOST_CHECK_EQUAL(prop1, 42);
      BOOST_CHECK_EQUAL(prop2, "3.141592");
      BOOST_CHECK_EQUAL(prop3, "nested_value");
    }
  };

  setenv("VOTCASHARE", ".", 1);
  char buff[FILENAME_MAX];
  std::cout << "WARNING: the VOTCASHARE env. variable has been updated to "
            << getcwd(buff, FILENAME_MAX) << "\n";

  // Generate user options
  tools::Property user_options;

  tools::Property &opt = user_options.add("options", "");
  tools::Property &opt_test = opt.add("testcalc", "");
  opt_test.add("option1", "42");
  tools::Property &new_prop = opt_test.add("option3", "");
  new_prop.add("nested", "nested_value");

  TestCalc test_calc;
  test_calc.Initialize(user_options);
}

BOOST_AUTO_TEST_SUITE_END()
