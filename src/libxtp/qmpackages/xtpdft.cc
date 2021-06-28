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
  log_file_name_ = job_name + ".orb";
  mo_file_name_ = log_file_name_;
  xtpdft_options_ = ParseCommonOptions(options);
}

bool XTPDFT::WriteInputFile(const Orbitals& orbitals) {
  orbitals_ = orbitals;
  orbitals_.setQMpackage(getPackageName());
  return true;
}

/**
 * Run calls DFTENGINE
 */
bool XTPDFT::RunDFT() {
  DFTEngine xtpdft;
  xtpdft.Initialize(xtpdft_options_);
  xtpdft.setLogger(pLog_);

  if (settings_.get<bool>("write_charges")) {
    xtpdft.setExternalcharges(&externalsites_);
  }
  bool success = xtpdft.Evaluate(orbitals_);
  std::string file_name = run_dir_ + "/" + log_file_name_;
  XTP_LOG(Log::error, *pLog_)
      << "Writing result to " << log_file_name_ << flush;
  orbitals_.WriteToCpt(file_name);
  return success;
}

void XTPDFT::CleanUp() {
  if (cleanup_.size() != 0) {
    XTP_LOG(Log::info, *pLog_) << "Removing " << cleanup_ << " files" << flush;
    std::vector<std::string> cleanup_info =
        tools::Tokenizer(cleanup_, ", ").ToVector();
    for (const std::string& substring : cleanup_info) {
      if (substring == "log") {
        std::string file_name = run_dir_ + "/" + log_file_name_;
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
    std::string file_name = run_dir_ + "/" + log_file_name_;
    orbitals.ReadFromCpt(file_name);
    XTP_LOG(Log::error, *pLog_) << (boost::format("QM energy[Hrt]: %4.8f ") %
                                    orbitals.getDFTTotalEnergy())
                                       .str()
                                << flush;
  } catch (std::runtime_error& error) {
    XTP_LOG(Log::error, *pLog_)
        << "Reading" << log_file_name_ << " failed" << flush;
    return false;
  }
  return true;
}

}  // namespace xtp
}  // namespace votca
