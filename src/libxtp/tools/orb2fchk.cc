/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

#include "orb2fchk.h"

// Third party includes
#include <boost/algorithm/string.hpp>

// Local VOTCA includes
#include "votca/xtp/gaussianwriter.h"
#include "votca/xtp/orbitals.h"
#include <votca/tools/constants.h>

namespace votca {
namespace xtp {

void Orb2Fchk::ParseOptions(const tools::Property& options) {

  basename_ = job_name_;
  orbfile_ = job_name_ + ".orb";
  state_ = options.get(".qmstate").as<QMState>();
  diff2gs_ = options.get(".diff2gs").as<bool>();
}

bool Orb2Fchk::Run() {
  log_.setReportLevel(Log::current_level);
  log_.setMultithreading(true);
  log_.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, log_) << "Loading data from " << orbfile_ << std::flush;
  XTP_LOG(Log::error, log_)
      << "Using density of state:  " << state_.ToString() << std::flush;
  orbitals.ReadFromCpt(orbfile_);

  GaussianWriter writer(log_);
  writer.WriteFile(basename_, orbitals, state_, diff2gs_);

  return true;
}

}  // namespace xtp
}  // namespace votca
