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

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/dftengine.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/qmpackagefactory.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/staticregion.h"

// Local private VOTCA includes
#include "qminqm.h"


namespace votca {
namespace xtp {

void QMinQM::ParseOptions(const tools::Property& options) {
  // lets get the archive file name from the xyz file name
  archive_file_ = job_name_ + ".orb";
  options_ = options;
}

bool QMinQM::Run() {

  std::cout << "Setup works" << std::endl;

  log_.setReportLevel(Log::current_level);

  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  // Get orbitals object
  Orbitals orbitals;
  orbitals.ReadFromCpt(archive_file_);
  DFTEngine xtpdft;
  xtpdft.Initialize(options_);
  xtpdft.setLogger(&log_);
  bool success = xtpdft.EvaluateActiveRegion(orbitals);
  orbitals.WriteToCpt("qminqm");
  return success;
}

}  // namespace xtp
}  // namespace votca
