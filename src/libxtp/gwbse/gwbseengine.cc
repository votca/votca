/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmpackage.h>

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

  _archive_file = archive_filename;
  std::string key = Identify();

  std::string tasks_string = options.get(".tasks").as<std::string>();

  if (tasks_string.find("guess") != std::string::npos) {
    _do_guess = true;
  }
  if (tasks_string.find("input") != std::string::npos) {
    _do_dft_input = true;
  }
  if (tasks_string.find("dft") != std::string::npos) {
    _do_dft_run = true;
  }
  if (tasks_string.find("parse") != std::string::npos) {
    _do_dft_parse = true;
  }
  if (tasks_string.find("gwbse") != std::string::npos) {
    _do_gwbse = true;
  }

  // XML option file for GWBSE
  if (_do_gwbse) {
    std::string _gwbse_xml = options.get(".gwbse_options").as<std::string>();
    _gwbse_options.LoadFromXML(_gwbse_xml);
  }
  // DFT log and MO file names
  _MO_file = options.get(".mofile").as<std::string>();
  _dftlog_file = options.get(".dftlog").as<std::string>();

  // Logger redirection
  _redirect_logger = options.ifExistsReturnElseReturnDefault<bool>(
      ".redirect_logger", _redirect_logger);
  _logger_file = "gwbse.log";

  // for requested merged guess, two archived orbitals objects are needed
  if (_do_guess) {
    _guess_archiveA =
        options.ifExistsReturnElseThrowRuntimeError<std::string>(".archiveA");
    _guess_archiveB =
        options.ifExistsReturnElseThrowRuntimeError<std::string>(".archiveB");
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
  Logger gwbse_engine_logger(_pLog->getReportLevel());
  Logger* logger = _pLog;
  if (_redirect_logger) {
    gwbse_engine_logger.setMultithreading(false);
    gwbse_engine_logger.setPreface(logINFO, (format("\n ...")).str());
    gwbse_engine_logger.setPreface(logERROR, (format("\n ...")).str());
    gwbse_engine_logger.setPreface(logWARNING, (format("\n ...")).str());
    gwbse_engine_logger.setPreface(logDEBUG, (format("\n ...")).str());
    logger = &gwbse_engine_logger;
  }
  _qmpackage->setLog(logger);
  if (_do_dft_input) {
    // required for merged guess
    if (_qmpackage->GuessRequested() && _do_guess) {  // do not want to do an
                                                      // SCF loop for a dimer
      XTP_LOG_SAVE(logINFO, *logger)
          << "Guess requested, reading molecular orbitals" << flush;
      Orbitals orbitalsA, orbitalsB;
      orbitalsA.ReadFromCpt(_guess_archiveA);
      orbitalsB.ReadFromCpt(_guess_archiveB);
      orbitals.PrepareDimerGuess(orbitalsA, orbitalsB);
    }
    _qmpackage->WriteInputFile(orbitals);
  }
  if (_do_dft_run) {
    bool run_success = _qmpackage->Run();
    if (!run_success) {
      throw std::runtime_error("\n DFT-run failed. Stopping!");
    }
  }

  // parse DFT data, if required
  if (_do_dft_parse) {
    XTP_LOG_SAVE(logINFO, *logger) << "Parsing DFT data from " << _dftlog_file
                                   << " and " << _MO_file << flush;
    _qmpackage->setLogFileName(_dftlog_file);
    _qmpackage->setMOsFileName(_MO_file);

    bool Logfile_parse = _qmpackage->ParseLogFile(orbitals);
    if (!Logfile_parse) {
      throw std::runtime_error("\n Parsing DFT logfile " + _dftlog_file +
                               " failed. Stopping!");
    }
    bool Orbfile_parse = _qmpackage->ParseMOsFile(orbitals);
    if (!Orbfile_parse) {
      throw std::runtime_error("\n Parsing DFT orbfile " + _MO_file +
                               " failed. Stopping!");
    }
    _qmpackage->CleanUp();
  }

  // if no parsing of DFT data is requested, reload serialized orbitals object
  if (!_do_dft_parse && _do_gwbse) {
    XTP_LOG_SAVE(logINFO, *logger)
        << "Loading serialized data from " << _archive_file << flush;
    orbitals.ReadFromCpt(_archive_file);
  }
  tools::Property& output_summary = _summary.add("output", "");
  if (_do_gwbse) {
    GWBSE gwbse = GWBSE(orbitals);
    gwbse.setLogger(logger);
    gwbse.Initialize(_gwbse_options);
    gwbse.Evaluate();
    gwbse.addoutput(output_summary);
  }
  if (_redirect_logger) {
    WriteLoggerToFile(logger);
  }
  return;
}

void GWBSEEngine::WriteLoggerToFile(Logger* pLog) {
  std::ofstream ofs;
  ofs.open(_logger_file, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("Bad file handle: " + _logger_file);
  }
  ofs << (*pLog) << std::endl;
  ofs.close();
  return;
}

}  // namespace xtp
}  // namespace votca
