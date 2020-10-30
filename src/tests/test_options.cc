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

#define BOOST_TEST_MODULE options_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/options.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(options_test)

BOOST_AUTO_TEST_CASE(resolving_links) {

  std::string default_loc = std::string(XTP_TEST_DATA_FOLDER) + "/options/";
  Options opt(default_loc);

  votca::tools::Property prop;
  prop.LoadFromXML(std::string(XTP_TEST_DATA_FOLDER) + "/options/option.xml");
  opt.ResolveLinks(prop);

  BOOST_CHECK_EQUAL(prop.get("calc.optb.package.distance.go").as<std::string>(),
                    "yes");
}

BOOST_AUTO_TEST_SUITE_END()
