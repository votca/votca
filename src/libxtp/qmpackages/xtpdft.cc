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
#include "xtpdft.h"

namespace votca {
namespace xtp {
using namespace std;

void XTPDFT::Initialize(const tools::Property& options) {
  const std::string job_name =
      options.ifExistsReturnElseReturnDefault<std::string>("job_name", "votca");
  _log_file_name = job_name + ".orb";
  _mo_file_name = _log_file_name;
  _xtpdft_options = ParseCommonOptions(options);
}

bool XTPDFT::WriteInputFile(const Orbitals& orbitals) {
  _orbitals = orbitals;
  _orbitals.setQMpackage(getPackageName());
  return true;
}

/**
 * Run calls DFTENGINE
 */
bool XTPDFT::Run() {
  DFTEngine xtpdft;
  xtpdft.Initialize(_xtpdft_options);
  xtpdft.setLogger(_pLog);

  if (_settings.get<bool>("write_charges")) {
    xtpdft.setExternalcharges(&_externalsites);
  }
  bool success = xtpdft.Evaluate(_orbitals);
  std::string file_name = _run_dir + "/" + _log_file_name;
  XTP_LOG(Log::error, *_pLog)
      << "Writing result to " << _log_file_name << flush;
  _orbitals.WriteToCpt(file_name);
  return success;
}

void XTPDFT::CleanUp() {
  if (_cleanup.size() != 0) {
    XTP_LOG(Log::info, *_pLog) << "Removing " << _cleanup << " files" << flush;
    tools::Tokenizer tok_cleanup(_cleanup, ", ");
    std::vector<std::string> cleanup_info;
    tok_cleanup.ToVector(cleanup_info);
    for (const std::string& substring : cleanup_info) {
      if (substring == "log") {
        std::string file_name = _run_dir + "/" + _log_file_name;
        remove(file_name.c_str());
      }
    }
  }

  return;
}

/**
 * Dummy, because XTPDFT adds info to orbitals directly
 */
bool XTPDFT::ParseMOsFile(Orbitals&) { return true; }

bool XTPDFT::ParseLogFile(Orbitals& orbitals) {
  try {
    std::string file_name = _run_dir + "/" + _log_file_name;
    orbitals.ReadFromCpt(file_name);
    XTP_LOG(Log::error, *_pLog) << (boost::format("QM energy[Hrt]: %4.8f ") %
                                    orbitals.getDFTTotalEnergy())
                                       .str()
                                << flush;
  } catch (std::runtime_error& error) {
    XTP_LOG(Log::error, *_pLog)
        << "Reading" << _log_file_name << " failed" << flush;
    return false;
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
