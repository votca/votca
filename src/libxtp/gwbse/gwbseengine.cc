/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations


#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/gwbse.h>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <votca/ctp/logger.h>



using boost::format;
using namespace boost::filesystem;
using std::flush;
namespace votca {
    namespace xtp {
      
        // +++++++++++++++++++++++++++++ //
        // GWBSEENGINE MEMBER FUNCTIONS  //
        // +++++++++++++++++++++++++++++ //

        void GWBSEEngine::Initialize(tools::Property& options, std::string archive_filename) {


            _archive_file = archive_filename;
            std::string key = Identify();

            // get the tasks
            std::string _tasks_string = options.get(".tasks").as<std::string> ();
            _do_guess = false;
            _do_dft_input = false;
            _do_dft_run = false;
            _do_dft_parse = false;
            _do_gwbse = false;

            if (_tasks_string.find("guess") != std::string::npos) _do_guess = true;
            if (_tasks_string.find("input") != std::string::npos) _do_dft_input = true;
            if (_tasks_string.find("dft") != std::string::npos)   _do_dft_run = true;
            if (_tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
            if (_tasks_string.find("gwbse") != std::string::npos) _do_gwbse = true;

            // XML option file for GWBSE
            std::string _gwbse_xml = options.get(".gwbse_options").as<std::string> ();
            load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());

            // DFT log and MO file names
            _MO_file = options.get(".mofile").as<std::string> ();
            _dftlog_file = options.get(".dftlog").as<std::string> ();

            // Logger redirection
            _redirect_logger = options.ifExistsReturnElseReturnDefault<bool>(".redirect_logger", false);
            _logger_file = "gwbse.log";
            
            // for requested merged guess, two archived orbitals objects are needed
            if ( _do_guess ){
                _guess_archiveA = options.ifExistsReturnElseThrowRuntimeError<std::string>(".archiveA");
                _guess_archiveB = options.ifExistsReturnElseThrowRuntimeError<std::string>(".archiveB");
            }

            return;
        }

        /* 
         *    CALL DFT and GWBSE modules to get excitation energies
         * 
         */


        void GWBSEEngine::ExcitationEnergies(QMPackage* qmpackage, Orbitals& orbitals) {

            //redirect log, if required
            // define own logger for GW-BSE that is written into a runFolder logfile
            ctp::Logger gwbse_engine_logger(_pLog->getReportLevel());
            ctp::Logger* logger=_pLog;
            if (_redirect_logger) {
                gwbse_engine_logger.setMultithreading(false);
                gwbse_engine_logger.setPreface(ctp::logINFO, (format("\n ...")).str());
                gwbse_engine_logger.setPreface(ctp::logERROR, (format("\n ...")).str());
                gwbse_engine_logger.setPreface(ctp::logWARNING, (format("\n ...")).str());
                gwbse_engine_logger.setPreface(ctp::logDEBUG, (format("\n ...")).str());
                logger=&gwbse_engine_logger;
            }
            qmpackage->setLog(logger);
            if (_do_dft_input) {
                // required for merged guess
                if (qmpackage->GuessRequested() && _do_guess) { // do not want to do an SCF loop for a dimer
                    CTP_LOG_SAVE(ctp::logINFO,*logger) << "Guess requested, reading molecular orbitals" << flush;
                    Orbitals orbitalsA, orbitalsB;
                    orbitalsA.ReadFromCpt(_guess_archiveA);
                    orbitalsB.ReadFromCpt(_guess_archiveB);
                    Orbitals::PrepareGuess(orbitalsA, orbitalsB, orbitals);
                }  
                qmpackage->WriteInputFile(orbitals);
            }
            if (_do_dft_run) {
                bool run_success = qmpackage->Run( orbitals );
                if (!run_success) {
                    throw std::runtime_error(std::string("\n DFT-run failed. Stopping!"));
                }
            }

            // parse DFT data, if required
            if (_do_dft_parse && qmpackage->getPackageName()!="xtp") {
                  
                     CTP_LOG_SAVE(ctp::logINFO,*logger) << "Parsing DFT data from " << _dftlog_file << " and " << _MO_file << flush;
                    qmpackage->setLogFileName(_dftlog_file);
                    qmpackage->setOrbitalsFileName(_MO_file);
                    qmpackage->ParseLogFile(orbitals);
                    qmpackage->ParseOrbitalsFile(orbitals);
            }

            // if no parsing of DFT data is requested, reload serialized orbitals object
            if (!_do_dft_parse && _do_gwbse) {
                CTP_LOG_SAVE(ctp::logINFO, *logger) << "Loading serialized data from " << _archive_file << flush;
                orbitals.ReadFromCpt(_archive_file);
            }

            if (_do_gwbse) {
                GWBSE gwbse = GWBSE(orbitals);
                gwbse.setLogger(logger);
                gwbse.Initialize(_gwbse_options);
                gwbse.Evaluate();
                if (_redirect_logger) WriteLoggerToFile(logger);
                tools::Property &output_summary = _summary.add("output", "");
                gwbse.addoutput(output_summary);
            }
            return;
        }

        void GWBSEEngine::WriteLoggerToFile(ctp::Logger* pLog) {
            std::ofstream ofs;
            ofs.open(_logger_file.c_str(), std::ofstream::out);
            if (!ofs.is_open()) {
                throw std::runtime_error("Bad file handle: " + _logger_file);
            }
            ofs << (*pLog) << std::endl;
            ofs.close();
            return;
        }

    }
}
