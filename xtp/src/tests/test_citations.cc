/*
 *            Copyright 2009-2026 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE citations_test

// Standard includes
#include <string>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/citations.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(citations_test)

BOOST_AUTO_TEST_CASE(table_keys_exist) {
  for (const std::string& name : citations::Calculators()) {
    const citations::CalculatorEntry entry = citations::ForCalculator(name);
    for (const std::string& key : entry.keys) {
      BOOST_CHECK_MESSAGE(citations::Exists(key), name + ": " + key);
    }
    for (const citations::Optional& opt : entry.optional) {
      BOOST_CHECK_MESSAGE(citations::Exists(opt.key), name + ": " + opt.key);
    }
  }
  BOOST_CHECK(citations::ForCalculator("nosuchcalculator").keys.empty());
}

BOOST_AUTO_TEST_CASE(unique_keys) {
  const std::vector<citations::Reference>& refs = citations::References();
  for (std::size_t i = 0; i < refs.size(); ++i) {
    for (std::size_t j = i + 1; j < refs.size(); ++j) {
      BOOST_CHECK(refs[i].key != refs[j].key);
    }
  }
}

BOOST_AUTO_TEST_CASE(cite_and_report) {
  citations::Reset();
  BOOST_CHECK(citations::Report().empty());
  citations::Cite("Baumeier:2024");
  citations::Cite("Chen:2025a");
  citations::Cite("Baumeier:2024");
  const std::vector<std::string> used = citations::Used();
  BOOST_REQUIRE_EQUAL(used.size(), 2);
  BOOST_CHECK_EQUAL(used[0], "Baumeier:2024");
  BOOST_CHECK_EQUAL(used[1], "Chen:2025a");
  const std::string report = citations::Report();
  BOOST_CHECK(report.find("[1]") != std::string::npos);
  BOOST_CHECK(report.find("doi:10.21105/joss.06864") != std::string::npos);
  BOOST_CHECK(report.find("[2]") != std::string::npos);
  BOOST_CHECK(report.find("[3]") == std::string::npos);
  citations::Reset();
  BOOST_CHECK(citations::Used().empty());
}

BOOST_AUTO_TEST_CASE(calculator_report) {
  const std::string text = citations::CalculatorReport("grw");
  BOOST_CHECK(text.find("10.1137/24M1658085") != std::string::npos);
  BOOST_CHECK(text.find("Depending") == std::string::npos);
  const std::string qmmm = citations::CalculatorReport("qmmm");
  BOOST_CHECK(qmmm.find("Depending on the options") != std::string::npos);
}

BOOST_AUTO_TEST_SUITE_END()
