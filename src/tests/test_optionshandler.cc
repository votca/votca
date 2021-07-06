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

#include <boost/test/tools/old/interface.hpp>
#include <stdexcept>
#define BOOST_TEST_MAIN

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/optionshandler.h"

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(optionshandler_test)

BOOST_AUTO_TEST_CASE(resolving_links) {

  std::string default_loc =
      std::string(TOOLS_TEST_DATA_FOLDER) + "/optionshandler/";

  std::string calcname = "testcalc";
  OptionsHandler opt(default_loc);
  Property result = opt.CalculatorOptions(calcname);

  BOOST_CHECK(result.exists("options.testcalc.subpackage.c"));
  BOOST_CHECK(result.exists("options.testcalc.subpackage.levA"));
  BOOST_CHECK(result.exists("options.testcalc.subpackage.levB"));
  BOOST_CHECK(!result.get("options.testcalc.subpackage").hasAttribute("link"));
  BOOST_CHECK(result.exists("options.testcalc.subpackage.subsubpackage.levB"));
}

BOOST_AUTO_TEST_CASE(check_required) {

  std::string default_loc =
      std::string(TOOLS_TEST_DATA_FOLDER) + "/optionshandler/";

  std::string calcname = "calc_required";
  OptionsHandler opt(default_loc);

  Property user_input;
  user_input.addTree("options." + calcname + ".a", "3");
  Property result = opt.ProcessUserInput(user_input, calcname);

  BOOST_CHECK(result.get("options.calc_required.a").as<votca::Index>() == 3);
  Property user_input2;
  user_input2.addTree("options." + calcname, "1");
  BOOST_CHECK_THROW(opt.ProcessUserInput(user_input2, calcname),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(check_optional) {

  std::string default_loc =
      std::string(TOOLS_TEST_DATA_FOLDER) + "/optionshandler/";

  std::string calcname = "calc_optional";
  OptionsHandler opt(default_loc);

  Property user_input;
  user_input.addTree("options." + calcname + ".a", "3");
  Property result = opt.ProcessUserInput(user_input, calcname);

  BOOST_CHECK(result.get("options.calc_optional.a").as<votca::Index>() == 3);
  Property user_input2;
  user_input2.addTree("options." + calcname, "1");
  BOOST_REQUIRE(!opt.ProcessUserInput(user_input2, calcname)
                     .exists("options.calc_optional.a"));
}

BOOST_AUTO_TEST_CASE(check_choices) {

  std::string default_loc =
      std::string(TOOLS_TEST_DATA_FOLDER) + "/optionshandler/";

  std::string calcname = "calc_choices";
  OptionsHandler opt(default_loc);

  Property user_input;
  user_input.addTree("options." + calcname, "");
  Property result = opt.ProcessUserInput(user_input, calcname);
  BOOST_CHECK(result.get("options.calc_choices.g").as<std::string>() == "a");

  Property user_input1;
  user_input1.addTree("options." + calcname + ".tasks", "input,dft,blabla");
  BOOST_CHECK_THROW(opt.ProcessUserInput(user_input1, calcname),
                      std::runtime_error);

  Property user_input_falseoption;
  user_input_falseoption.addTree("options." + calcname + ".oinkoink", "input,dft,blabla");
  BOOST_CHECK_THROW(opt.ProcessUserInput(user_input_falseoption, calcname),
                      std::runtime_error);

  Property user_input2;
  user_input2.addTree("options." + calcname + ".g", "e");
  BOOST_CHECK_THROW(opt.ProcessUserInput(user_input2, calcname),
                      std::runtime_error);

  Property user_input3;
  user_input3.addTree("options." + calcname + ".b", "-1");
  BOOST_CHECK_THROW(opt.ProcessUserInput(user_input3, calcname),
                      std::runtime_error);

  Property user_input4;
  user_input4.addTree("options." + calcname + ".c", "a");
  BOOST_CHECK_THROW(opt.ProcessUserInput(user_input4, calcname),
                      std::runtime_error);

  Property user_input5;
  user_input5.addTree("options." + calcname + ".d", "a");
  BOOST_CHECK_THROW(opt.ProcessUserInput(user_input5, calcname),
                      std::runtime_error);
}


BOOST_AUTO_TEST_CASE(check_list) {

  std::string default_loc =
      std::string(TOOLS_TEST_DATA_FOLDER) + "/optionshandler/";

  std::string calcname = "calc_list";
  OptionsHandler opt(default_loc);

  Property user_input;
  user_input.addTree("options." + calcname + ".d.a", "3");
  user_input.addTree("options." + calcname + ".d.a", "4");
  user_input.addTree("options." + calcname + ".d.b.type", "4");
  auto& temp=user_input.addTree("options." + calcname + ".d.b","");
  temp.add("type","5");
  Property result = opt.ProcessUserInput(user_input, calcname);

  BOOST_CHECK(result.Select("options.calc_list.d.a").size() == 2);
  BOOST_CHECK(result.Select("options.calc_list.d.b").size() == 2);
  BOOST_CHECK(result.Select("options.calc_list.d.f").size() == 0);
  BOOST_CHECK(result.Select("options.calc_list.c").size() == 1);
}

BOOST_AUTO_TEST_CASE(check_brokenlist) {

  std::string default_loc =
      std::string(TOOLS_TEST_DATA_FOLDER) + "/optionshandler/";

  std::string calcname = "calc_brokenlist";
  OptionsHandler opt(default_loc);

  Property user_input;
  user_input.addTree("options." + calcname + ".d.a", "3");
  BOOST_CHECK_THROW(opt.ProcessUserInput(user_input, calcname),std::runtime_error);

}

BOOST_AUTO_TEST_SUITE_END()
