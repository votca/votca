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

#ifndef _VOTCA_XTP_GWBSEENGINE_H
#define _VOTCA_XTP_GWBSEENGINE_H

#include <votca/xtp/segment.h>
#include <votca/xtp/polarseg.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/apolarsite.h>
#include <boost/filesystem.hpp>
#include <votca/xtp/logger.h>

namespace votca {
    namespace xtp {
        class QMPackage;
        class Orbitals;

/**
         * \brief Electronic Excitations via Density-Functional Theory
         *
         * Evaluates electronic ground state in molecular systems based on
         * density functional theory with Gaussian Orbitals.
         * 
         */

        class GWBSEEngine {
        public:

            std::string Identify() {
                return "gwbse_engine";
            }

            void Initialize(tools::Property &options, std::string archive_filename);
            void ExcitationEnergies(Orbitals& orbitals);

            void setLog(xtp::Logger* pLog) {
                _pLog = pLog;
            }
            
            void setQMPackage(QMPackage* qmpackage){
                _qmpackage=qmpackage;
            }

            std::string GetDFTLog() {
                return _dftlog_file;
            };

            void setLoggerFile(std::string logger_file) {
                _logger_file = logger_file;
            };

            void setRedirectLogger(bool redirect_logger) {
                _redirect_logger = redirect_logger;
            };
            
            
            tools::Property& ReportSummary(){ return _summary;};


        private:
            
            QMPackage* _qmpackage;

            xtp::Logger *_pLog;

            // task options
            bool _do_guess;
            bool _do_dft_input;
            bool _do_dft_run;
            bool _do_dft_parse;
            bool _do_gwbse;
            bool _redirect_logger;

            // DFT log and MO file names
            std::string _MO_file; // file containing the MOs from qmpackage...
            std::string _dftlog_file; // file containing the Energies etc... from qmpackage...
            std::string _logger_file;
            std::string _archive_file;
            std::string _guess_archiveA;
            std::string _guess_archiveB;

            // Options for GWBSE module
            tools::Property _gwbse_options;
            tools::Property _summary;

            void WriteLoggerToFile(xtp::Logger* pLog);
        };


    }
}

#endif /* _VOTCA_XTP_GWBSEENGINE_H */
