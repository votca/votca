/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE dfttimings_test

// Standard includes
#include <sstream>
#include <thread>
#include <utility>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/dfttimings.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(dfttimings_test)

namespace {
void Sleep(double seconds) {
  std::this_thread::sleep_for(std::chrono::duration<double>(seconds));
}

// Reads the calls and time columns of one row of the report.
std::pair<long, double> ReportedRow(const std::string& report,
                                    const std::string& name) {
  std::size_t pos = report.find(name);
  BOOST_REQUIRE(pos != std::string::npos);
  std::istringstream line(report.substr(pos + name.size()));
  long calls = 0;
  double seconds = 0.0;
  line >> calls >> seconds;
  return {calls, seconds};
}
}  // namespace

// Each entry holds its own time only: nested scopes and added time are taken
// out of the enclosing scope, so nothing is counted twice.
BOOST_AUTO_TEST_CASE(nested_scopes_are_counted_once) {
  DFTTimings timings;
  timings.Reset();
  {
    auto outer = timings.Measure("outer");
    Sleep(0.05);
    {
      auto inner = timings.Measure("inner");
      Sleep(0.10);
    }
    {
      auto inner = timings.Measure("inner");
      Sleep(0.10);
    }
    timings.Add("added", 0.03);
  }

  Logger log;
  log.setReportLevel(votca::Log::error);
  timings.Report(log, votca::Log::error);
  std::stringstream ss;
  ss << log;
  const std::string report = ss.str();

  // outer: own 0.05 s, minus the 0.03 s added inside it
  BOOST_CHECK_EQUAL(ReportedRow(report, "inner").first, 2);
  BOOST_CHECK_CLOSE(ReportedRow(report, "inner").second, 0.20, 15);
  BOOST_CHECK_SMALL(ReportedRow(report, "outer").second - 0.02, 0.02);
  BOOST_CHECK_CLOSE(ReportedRow(report, "added").second, 0.03, 1);
}

BOOST_AUTO_TEST_CASE(resident_memory_is_read) {
  double rss = DFTTimings::ResidentMemoryGB(false);
  double peak = DFTTimings::ResidentMemoryGB(true);
#ifdef __linux__
  BOOST_CHECK_GT(rss, 0.0);
  BOOST_CHECK_GE(peak, rss);
#else
  BOOST_CHECK(rss < 0 || peak >= rss);
#endif
}

BOOST_AUTO_TEST_SUITE_END()
