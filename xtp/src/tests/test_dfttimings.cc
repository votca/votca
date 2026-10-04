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
#include <chrono>
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
// out of the enclosing scope, so nothing is counted twice. The expected
// values are measured around the scopes rather than taken from the sleep
// durations: sleep_for may oversleep a lot on a loaded machine (CI runners).
BOOST_AUTO_TEST_CASE(nested_scopes_are_counted_once) {
  using Clock = std::chrono::steady_clock;
  auto Seconds = [](Clock::time_point from, Clock::time_point to) {
    return std::chrono::duration<double>(to - from).count();
  };
  DFTTimings timings;
  timings.Reset();
  double inner_measured = 0.0;
  double outer_measured = 0.0;
  {
    const Clock::time_point outer_start = Clock::now();
    {
      auto outer = timings.Measure("outer");
      Sleep(0.05);
      for (int i = 0; i < 2; ++i) {
        const Clock::time_point start = Clock::now();
        {
          auto inner = timings.Measure("inner");
          Sleep(0.10);
        }
        inner_measured += Seconds(start, Clock::now());
      }
      timings.Add("added", 0.03);
    }
    outer_measured = Seconds(outer_start, Clock::now());
  }

  Logger log;
  log.setReportLevel(votca::Log::error);
  timings.Report(log, votca::Log::error);
  std::stringstream ss;
  ss << log;
  const std::string report = ss.str();

  // The report prints 2 decimals; the measured intervals enclose the scopes,
  // so they can only be slightly longer than what the scopes record.
  const double tolerance = 0.006 + 0.01;
  BOOST_CHECK_EQUAL(ReportedRow(report, "inner").first, 2);
  const double inner = ReportedRow(report, "inner").second;
  BOOST_CHECK_GE(inner, 0.20 - 0.006);
  BOOST_CHECK_SMALL(inner - inner_measured, tolerance);
  // outer: its own time, without the inner scopes and the 0.03 s added
  const double outer = ReportedRow(report, "outer").second;
  BOOST_CHECK_SMALL(outer - (outer_measured - inner_measured - 0.03),
                    tolerance);
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
