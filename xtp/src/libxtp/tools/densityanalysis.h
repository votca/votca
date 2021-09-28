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
#ifndef VOTCA_XTP_DENSITYANALYSIS_H
#define VOTCA_XTP_DENSITYANALYSIS_H

// Standard includes
#include <cstdio>

// Third party includes
#include <boost/filesystem.hpp>

// Local VOTCA includes
#include "votca/xtp/gyration.h"
#include "votca/xtp/logger.h"

namespace votca {
namespace xtp {

class DensityAnalysis final : public QMTool {
 public:
  std::string Identify() const { return "densityanalysis"; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Run();

 private:
  std::string orbfile_;
  std::string output_file_;
  tools::Property gyration_options_;

  Logger log_;
};

void DensityAnalysis::ParseOptions(const tools::Property& options) {

  orbfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".input", job_name_ + ".orb");

  gyration_options_ = options;
}

bool DensityAnalysis::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);

  log_.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, log_)
      << " Loading QM data from " << orbfile_ << std::flush;
  orbitals.ReadFromCpt(orbfile_);

  Density2Gyration density2gyration(log_);
  density2gyration.Initialize(gyration_options_);
  density2gyration.AnalyzeDensity(orbitals);

  return true;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DENSITYANALYSIS_H
