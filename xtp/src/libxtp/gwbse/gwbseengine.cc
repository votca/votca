/*
 *            Copyright 2009-2022 The VOTCA Development Team
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
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <string>

// Local VOTCA includes
#include "votca/xtp/gwbse.h"
#include "votca/xtp/gwbseengine.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/pmlocalization.h"
#include "votca/xtp/qmpackage.h"

using boost::format;
using namespace boost::filesystem;
using std::flush;
namespace votca {
namespace xtp {

// +++++++++++++++++++++++++++++ //
// GWBSEENGINE MEMBER FUNCTIONS  //
// +++++++++++++++++++++++++++++ //

void GWBSEEngine::Initialize(tools::Property& options,
                             std::string archive_filename) {

  archive_file_ = archive_filename;

  std::string tasks_string = options.get(".tasks").as<std::string>();

  // We split either on a space or a comma
  tools::Tokenizer tokenizedTasks(tasks_string, " ,");
  std::vector<std::string> tasks = tokenizedTasks.ToVector();
  // Check which tasks exist and set them to true
  do_guess_ = std::find(tasks.begin(), tasks.end(), "guess") != tasks.end();
  do_dft_input_ = std::find(tasks.begin(), tasks.end(), "input") != tasks.end();
  do_dft_run_ = std::find(tasks.begin(), tasks.end(), "dft") != tasks.end();
  do_dft_parse_ = std::find(tasks.begin(), tasks.end(), "parse") != tasks.end();
  do_gwbse_ = std::find(tasks.begin(), tasks.end(), "gwbse") != tasks.end();
  do_localize_ =
      std::find(tasks.begin(), tasks.end(), "localize") != tasks.end();
  do_dft_in_dft_ =
      std::find(tasks.begin(), tasks.end(), "dft_in_dft") != tasks.end();

  // XML option file for GWBSE
  if (do_gwbse_) {
    gwbse_options_ = options.get(".gwbse");
  }
  if (do_localize_) {
    localize_options_ = options.get(".localize");
  }
  // DFT log and MO file names
  MO_file_ = qmpackage_->getMOFile();
  dftlog_file_ = qmpackage_->getLogFile();

  // Logger redirection
  if (options.exists(".logging_file")) {
    logger_file_ = options.get(".logging_file").as<std::string>();
  }

  // for requested merged guess, two archived orbitals objects are needed
  if (do_guess_) {
    guess_archiveA_ = options.get(".archiveA").as<std::string>();
    guess_archiveB_ = options.get(".archiveB").as<std::string>();
  }

  return;
}

/*
 *    CALL DFT and GWBSE modules to get excitation energies
 *
 */

void GWBSEEngine::ExcitationEnergies(Orbitals& orbitals) {

  // redirect log, if required
  // define own logger for GW-BSE that is written into a runFolder logfile
  Logger gwbse_engine_logger(pLog_->getReportLevel());
  Logger* logger = pLog_;
  if (!logger_file_.empty()) {
    gwbse_engine_logger.setMultithreading(false);
    gwbse_engine_logger.setPreface(Log::info, "\n ...");
    gwbse_engine_logger.setPreface(Log::error, "\n ...");
    gwbse_engine_logger.setPreface(Log::warning, "\n ...");
    gwbse_engine_logger.setPreface(Log::debug, "\n ...");
    logger = &gwbse_engine_logger;
  }
  qmpackage_->setLog(logger);
  if (do_dft_input_) {
    // required for merged guess
    if (qmpackage_->GuessRequested() && do_guess_) {  // do not want to do an
                                                      // SCF loop for a dimer
      XTP_LOG(Log::error, *logger)
          << "Guess requested, reading molecular orbitals" << flush;
      Orbitals orbitalsA, orbitalsB;
      orbitalsA.ReadFromCpt(guess_archiveA_);
      orbitalsB.ReadFromCpt(guess_archiveB_);
      orbitals.PrepareDimerGuess(orbitalsA, orbitalsB);
    }
    qmpackage_->WriteInputFile(orbitals);
  }
  if (do_dft_run_) {
    bool run_success = qmpackage_->Run();
    if (!run_success) {
      throw std::runtime_error("\n DFT-run failed. Stopping!");
    }
  }

  // parse DFT data, if required
  if (do_dft_parse_) {
    XTP_LOG(Log::error, *logger) << "Parsing DFT data from " << dftlog_file_
                                 << " and " << MO_file_ << flush;
    qmpackage_->setLogFileName(dftlog_file_);
    qmpackage_->setMOsFileName(MO_file_);

    bool Logfile_parse = qmpackage_->ParseLogFile(orbitals);
    if (!Logfile_parse) {
      throw std::runtime_error("\n Parsing DFT logfile " + dftlog_file_ +
                               " failed. Stopping!");
    }
    bool Orbfile_parse = qmpackage_->ParseMOsFile(orbitals);
    if (!Orbfile_parse) {
      throw std::runtime_error("\n Parsing DFT orbfile " + MO_file_ +
                               " failed. Stopping!");
    }
    qmpackage_->CleanUp();
  }

  if (do_localize_) {
    if (!do_dft_parse_) {
      orbitals.ReadFromCpt(archive_file_);
    }
    PMLocalization pml(*logger, localize_options_);
    pml.computePML(orbitals);
  }

  if (!do_dft_parse_ && do_gwbse_ || !do_dft_parse_ && do_dft_in_dft_) {
    XTP_LOG(Log::error, *logger)
        << "Loading serialized data from " << archive_file_ << flush;
    orbitals.ReadFromCpt(archive_file_);
  }

  if (do_dft_in_dft_ && !do_localize_) {
    if (orbitals.getLMOs().size() == 0) {
      throw std::runtime_error(
          "Can't do DFT in DFT embedding without localization");
    }
  }

  if (do_dft_in_dft_) {
    qmpackage_->WriteInputFile(orbitals);
    bool run_success = qmpackage_->RunActiveRegion();
    if (!run_success) {
      throw std::runtime_error("\n DFT in DFT embedding failed. Stopping!");
    }
    bool Logfile_parse = qmpackage_->ParseLogFile(orbitals);
    if (!Logfile_parse) {
      throw std::runtime_error("\n Parsing DFT logfile " + dftlog_file_ +
                               " failed. Stopping!");
    }
    bool Orbfile_parse = qmpackage_->ParseMOsFile(orbitals);
    if (!Orbfile_parse) {
      throw std::runtime_error("\n Parsing DFT orbfile " + MO_file_ +
                               " failed. Stopping!");
    }
  }

  // if no parsing of DFT data is requested, reload serialized orbitals object

  tools::Property& output_summary = summary_.add("output", "");

  if (do_gwbse_) {
    if (do_dft_in_dft_ || orbitals.getEmbeddedMOs().eigenvectors().size() > 0) {
      Orbitals orb_embedded = orbitals;
      orb_embedded.MOs() = orb_embedded.getEmbeddedMOs();
      Index active_electrons = orb_embedded.getNumOfActiveElectrons();
      orb_embedded.setNumberOfAlphaElectrons(active_electrons);
      orb_embedded.setNumberOfOccupiedLevels(active_electrons / 2);
      GWBSE gwbse = GWBSE(orb_embedded);
      gwbse.setLogger(logger);
      gwbse.Initialize(gwbse_options_);
      gwbse.Evaluate();
      gwbse.addoutput(output_summary);
      orbitals = orb_embedded;
    } else {
      GWBSE gwbse = GWBSE(orbitals);
      gwbse.setLogger(logger);
      gwbse.Initialize(gwbse_options_);
      gwbse.Evaluate();
      gwbse.addoutput(output_summary);
    }
  }

  if (!logger_file_.empty()) {
    WriteLoggerToFile(logger);
  }
  return;
}

void GWBSEEngine::WriteLoggerToFile(Logger* pLog) {
  std::ofstream ofs;
  ofs.open(logger_file_, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("Bad file handle: " + logger_file_);
  }
  ofs << (*pLog) << std::endl;
  ofs.close();
  return;
}

}  // namespace xtp
}  // namespace votca
