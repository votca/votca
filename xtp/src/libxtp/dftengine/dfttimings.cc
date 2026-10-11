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

// Third party includes
#include <boost/format.hpp>

// Local VOTCA includes
#include "votca/xtp/dfttimings.h"

namespace votca {
namespace xtp {

void DFTTimings::Report(Logger& log, Log::Level level,
                        const std::string& title) const {
  const double total = Elapsed();
  double accounted = 0.0;
  XTP_LOG(level, log) << " " << title << " (wall clock)" << std::flush;
  XTP_LOG(level, log) << (boost::format("   %-38s %7s %11s %10s %6s") % "part" %
                          "calls" % "time [s]" % "per call" % "%")
                             .str()
                      << std::flush;
  for (const Entry& e : entries_) {
    accounted += e.seconds;
    XTP_LOG(level, log) << (boost::format("   %-38s %7d %11.2f %10.3f %6.1f") %
                            e.name % e.calls % e.seconds %
                            (e.seconds / double(e.calls)) %
                            (total > 0 ? 100.0 * e.seconds / total : 0.0))
                               .str()
                        << std::flush;
  }
  XTP_LOG(level, log) << (boost::format("   %-38s %7s %11.2f %10s %6.1f") %
                          "not itemised" % "" % (total - accounted) % "" %
                          (total > 0 ? 100.0 * (total - accounted) / total
                                     : 0.0))
                             .str()
                      << std::flush;
  XTP_LOG(level, log)
      << (boost::format("   %-38s %7s %11.2f") % "total" % "" % total).str()
      << std::flush;
}

std::string DFTTimings::Format() const {
  std::string line;
  for (const Entry& e : entries_) {
    if (!line.empty()) {
      line += ", ";
    }
    line += (boost::format("%s %.2f s") % e.name % e.seconds).str();
  }
  return line;
}

double DFTTimings::ResidentMemoryGB(bool peak) {
  std::ifstream status("/proc/self/status");
  if (!status.is_open()) {
    return -1.0;
  }
  const std::string key = peak ? "VmHWM:" : "VmRSS:";
  std::string word;
  while (status >> word) {
    if (word == key) {
      double kb = 0.0;
      status >> kb;
      return kb / (1024.0 * 1024.0);
    }
  }
  return -1.0;
}

}  // namespace xtp
}  // namespace votca
