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

// Standard includes
#include <cstdio>
#include <iomanip>

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local private VOTCA includes
#include "votca/xtp/xtpdft.h"

namespace votca {
namespace xtp {
using namespace std;

void XTPDFT::ParseSpecificOptions(const tools::Property& options) {
  const std::string job_name = options.get("temporary_file").as<std::string>();
  log_file_name_ = job_name + ".orb";
  mo_file_name_ = log_file_name_;
}

bool XTPDFT::WriteInputFile(const Orbitals&) {
  return true; // dummy, no input files to write for xtp dft
}

/**
 * Run calls DFTENGINE
 */
bool XTPDFT::RunDFT(Orbitals& orbitals) {
  orbitals.setQMpackage(getPackageName());
  DFTEngine xtpdft;
  xtpdft.Initialize(options_);
  xtpdft.setLogger(pLog_);

  if (!externalsites_.empty()) {
    xtpdft.setExternalcharges(&externalsites_);
  }
  bool success = xtpdft.Evaluate(orbitals);
  if (!success) {
    return success;
  }

  XTP_LOG(Log::error, *pLog_) << (boost::format("QM energy[Hrt]: %4.8f ") %
                                  orbitals.getDFTTotalEnergy())
                                     .str()
                              << flush;
  return true;
}

}  // namespace xtp
}  // namespace votca
