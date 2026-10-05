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

#pragma once
#ifndef VOTCA_XTP_DFTTIMINGS_H
#define VOTCA_XTP_DFTTIMINGS_H

// Standard includes
#include <chrono>
#include <fstream>
#include <string>
#include <vector>

// Local VOTCA includes
#include "logger.h"

namespace votca {
namespace xtp {

/**
 * \brief Wall-clock accounting for the parts of a DFT run.
 *
 * Entries keep the order in which they were first used, so the report reads
 * like the run. Scopes may nest; each entry holds only its own time, without
 * the time of scopes nested inside it, so the entries add up to the total.
 * Scopes must be opened and closed on one thread (the serial driver code),
 * not inside OpenMP regions.
 */
class DFTTimings {
 public:
  using Clock = std::chrono::steady_clock;

  class Scope {
   public:
    Scope(DFTTimings& timings, std::string name)
        : timings_(timings), name_(std::move(name)), start_(Clock::now()) {
      timings_.nested_.push_back(0.0);
    }
    ~Scope() {
      double elapsed =
          std::chrono::duration<double>(Clock::now() - start_).count();
      double nested = timings_.nested_.back();
      timings_.nested_.pop_back();
      timings_.AddEntry(name_, elapsed - nested);
      if (!timings_.nested_.empty()) {
        timings_.nested_.back() += elapsed;
      }
    }
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;

   private:
    DFTTimings& timings_;
    std::string name_;
    Clock::time_point start_;
  };

  Scope Measure(const std::string& name) { return Scope(*this, name); }

  /// Adds time measured elsewhere. Inside an open scope it is taken out of
  /// that scope's own time, as if it had been a nested scope.
  void Add(const std::string& name, double seconds) {
    AddEntry(name, seconds);
    if (!nested_.empty()) {
      nested_.back() += seconds;
    }
  }

  void Reset() {
    entries_.clear();
    nested_.clear();
    start_ = Clock::now();
  }

  double Elapsed() const {
    return std::chrono::duration<double>(Clock::now() - start_).count();
  }

  /// Table of all entries with their share of the time since Reset().
  void Report(Logger& log, Log::Level level,
              const std::string& title = "DFT timing summary") const;

  /// Resident memory of this process in GB, read from /proc/self/status
  /// (VmRSS, or VmHWM for the peak). Returns a negative value where that
  /// file does not exist.
  static double ResidentMemoryGB(bool peak);

 private:
  void AddEntry(const std::string& name, double seconds) {
    for (Entry& e : entries_) {
      if (e.name == name) {
        e.seconds += seconds;
        ++e.calls;
        return;
      }
    }
    entries_.push_back({name, seconds, 1});
  }

  struct Entry {
    std::string name;
    double seconds;
    long calls;
  };
  std::vector<Entry> entries_;
  std::vector<double> nested_;  // time of finished inner scopes, per level
  Clock::time_point start_ = Clock::now();
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTTIMINGS_H
