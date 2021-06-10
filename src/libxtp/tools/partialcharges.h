/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
#ifndef VOTCA_XTP_PARTIALCHARGES_H
#define VOTCA_XTP_PARTIALCHARGES_H

// Standard includes
#include <cstdio>

// Third party includes
#include <boost/filesystem.hpp>

// Local VOTCA includes
#include "votca/xtp/esp2multipole.h"
#include "votca/xtp/logger.h"

namespace votca {
namespace xtp {

class Partialcharges final : public QMTool {
 public:
  Partialcharges() = default;
  ~Partialcharges() = default;

  std::string Identify() { return "partialcharges"; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Run();

 private:
  std::string orbfile_;
  std::string output_file_;
  tools::Property esp_options_;

  Logger log_;
};

void Partialcharges::ParseOptions(const tools::Property& options) {

  orbfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".input", job_name_ + ".orb");
  output_file_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".output", job_name_ + ".mps");
  esp_options_ = options.get(".esp_options");
}

bool Partialcharges::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);

  log_.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, log_)
      << " Loading QM data from " << orbfile_ << std::flush;
  orbitals.ReadFromCpt(orbfile_);
  Esp2multipole esp2multipole = Esp2multipole(log_);
  esp2multipole.Initialize(esp_options_);
  StaticSegment seg = esp2multipole.Extractingcharges(orbitals);
  seg.WriteMPS(output_file_, esp2multipole.GetStateString());

  XTP_LOG(Log::error, log_)
      << "Written charges to " << output_file_ << std::flush;

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_PARTIALCHARGES_H
