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

  _basename = _job_name;
  _orbfile = _job_name + ".orb";
  _state_string = options.get(".qmstate").as<std::string>();
}

bool Orb2Fchk::Run() {
  _log.setReportLevel(Log::current_level);
  _log.setMultithreading(true);
  _log.setCommonPreface("\n... ...");

  Orbitals orbitals;
  XTP_LOG(Log::error, _log) << "Loading data from " << _orbfile << std::flush;
  XTP_LOG(Log::error, _log)
      << "Using density of state:  " << _state_string << std::flush;
  orbitals.ReadFromCpt(_orbfile);

  GaussianWriter writer(_log);
  writer.WriteFile(_basename, orbitals, QMState(_state_string));

  return true;
}

}  // namespace xtp
}  // namespace votca
