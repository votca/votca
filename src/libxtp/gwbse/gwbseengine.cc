/*
 *            Copyright 2009-2017 The VOTCA Development Team
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
#include <boost/numeric/ublas/operation.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>

#include <votca/ctp/logger.h>



using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;
        // +++++++++++++++++++++++++++++ //
        // GWBSEENGINE MEMBER FUNCTIONS  //
        // +++++++++++++++++++++++++++++ //

        void GWBSEENGINE::Initialize(Property* options, string _archive_filename) {

            
            _archive_file = _archive_filename;
            string key = Identify();

            // get the tasks
            string _tasks_string = options->get(".tasks").as<string> ();
            _do_dft_input = false;
            _do_dft_run   = false;
            _do_dft_parse = false;
            _do_gwbse     = false;
            if (_tasks_string.find("input") != std::string::npos) _do_dft_input = true;
            if (_tasks_string.find("dft") != std::string::npos)   _do_dft_run   = true;
            if (_tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
            if (_tasks_string.find("gwbse") != std::string::npos) _do_gwbse     = true;
            
            // XML option file for GWBSE
            string _gwbse_xml = options->get(".gwbse_options").as<string> ();
            load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());

            // DFT log and MO file names
            _MO_file     = options->get(".mofile").as<string> ();
            _dftlog_file = options->get(".dftlog").as<string> ();
            
            
            return;
        }

        /* 
         *    CALL DFT and GWBSE modules to get excitation energies
         * 
         */


        void GWBSEENGINE::ExcitationEnergies(QMPackage* _qmpackage, vector<ctp::Segment*> _segments, Orbitals* _orbitals) {


            if (_do_dft_input) {
                _qmpackage->WriteInputFile(_segments);
            }

            if (_do_dft_run) {
                bool run_success = _qmpackage->Run();
                if (!run_success) {
                    throw runtime_error(string("\n GW-BSE without DFT is difficult. Stopping!"));
                }
            }

            // parse DFT data, if required
            if (_do_dft_parse) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Parsing DFT data from " << _dftlog_file << " and " << _MO_file << flush;
                _qmpackage->setLogFileName(_dftlog_file);
                _qmpackage->ParseLogFile(_orbitals);
                _qmpackage->setOrbitalsFileName(_MO_file);
                _qmpackage->ParseOrbitalsFile(_orbitals); 
                _orbitals->setDFTbasis(_qmpackage->getBasisSetName());
            }

            // if no parsing of DFT data is requested, reload serialized orbitals object
            if (!_do_dft_parse) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Loading serialized data from " << _archive_file << flush;
                _orbitals->Load(_archive_file);
            }

            if (_do_gwbse) {
                GWBSE _gwbse = GWBSE(_orbitals);
                _gwbse.setLogger(_pLog);
                _gwbse.Initialize(&_gwbse_options);
                _gwbse.Evaluate();
                
                // not sure what this does at the moment
                //Property *_output_summary = &(_summary.add("output", ""));
                //_gwbse.addoutput(_output_summary);
            }
            return;
        }


    }
}
