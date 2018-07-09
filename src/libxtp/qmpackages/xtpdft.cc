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

#include "xtpdft.h"
#include <votca/ctp/segment.h>
#include <votca/xtp/qminterface.h>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/constants.h>
#include <stdio.h>
#include <iomanip>
#include <sys/stat.h>
#include <vector>



namespace votca {
    namespace xtp {
      using namespace std;

        void XTPDFT::Initialize(tools::Property *options) {

            // GAUSSIAN file names
            std::string fileName = "system";

            _xyz_file_name = fileName + ".xyz";

            std::string key = "package";
            std::string _name = options->get(key + ".name").as<std::string> ();

            if (_name != "xtp") {
                cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error("Wrong options file");
            }

            _charge = options->get(key + ".charge").as<int> ();
            _spin = options->get(key + ".spin").as<int> ();
            _options = options->get(key + ".options").as<std::string> ();
            _threads = options->get(key + ".threads").as<int> ();
            _cleanup = options->get(key + ".cleanup").as<std::string> ();

            // pass the information about the dftengine options and init
            load_property_from_xml(_xtpdft_options, _options.c_str());
            _xtpdft.Initialize( &_xtpdft_options );

            // check if ECPs are used in xtpdft
            _write_pseudopotentials=false;
            if (_xtpdft_options.exists("dftengine.ecp")){
                if (_xtpdft_options.get("dftengine.ecp").as<std::string> () !="") {
                    _write_pseudopotentials=true;
                }
            }

        }

        /**
         * Dummy for use of XTPDFT as QMPackage, needs no input file
         */
        bool XTPDFT::WriteInputFile(std::vector<ctp::Segment* > segments, Orbitals* orbitals_guess, std::vector<ctp::PolarSeg* > polar_segments ) {

            //if ( orbitals_guess != NULL )  {
            //    CTP_LOG(ctp::logDEBUG, _log) << "Reading guess from " << _guess_file << flush;
            //    _orbitals.Load(_guess_file);

            //} else {

           // }
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Preparing XTP DFTENGINE "  << flush;
            _xtpdft.setLogger(_pLog);
            _xtpdft.Prepare( orbitals_guess );
            CTP_LOG(ctp::logDEBUG, *_pLog) << " done "  << flush;

            return true;
        }


        bool XTPDFT::setMultipoleBackground( std::vector<ctp::PolarSeg*> multipoles){

            _xtpdft.setExternalcharges(multipoles);

            return true;
        }



        /**
         * Run calls DFTENGINE
         */
        bool XTPDFT::Run( Orbitals* _orbitals ) {

            if ( !_orbitals->hasQMAtoms() ){
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Reading structure from " << _xyz_file_name << flush;
                _orbitals->LoadFromXYZ(_xyz_file_name);
            }


            CTP_LOG(ctp::logDEBUG, *_pLog) << "Running XTP DFT " << flush;
            //

            _xtpdft.Evaluate( _orbitals );
            _basisset_name = _xtpdft.getDFTBasisName();

            return true;

        }

        /**
         * Clean up dummy, may be required if use of scratch will be added
         */
        void XTPDFT::CleanUp() {

            return;
        }

        /**
         * Dummy, because XTPDFT adds info to orbitals directly
         */
        bool XTPDFT::ParseOrbitalsFile(Orbitals * _orbitals) {
            return true;
        }

        /* DON'T THINK THIS IS NEEDED */
        bool XTPDFT::CheckLogFile() {

            return true;
        }

        /**
         * Dummy, because information is directly stored in orbitals
         */
        bool XTPDFT::ParseLogFile(Orbitals * _orbitals) {

            return true;
        }



    }
}
