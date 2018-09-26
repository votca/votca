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
#include <votca/xtp/segment.h>
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

        void XTPDFT::Initialize(tools::Property &options) {
            _xtpdft_options=options;
            _log_file_name="system.orb";
            std::string key = "package";
            std::string packagename = _xtpdft_options.get(key + ".name").as<std::string> ();

            if (packagename != "xtp") {
                cerr << "Tried to use " << packagename << " package. ";
                throw std::runtime_error("Wrong options file");
            }

            _charge = _xtpdft_options.get(key + ".charge").as<int> ();
            _spin = _xtpdft_options.get(key + ".spin").as<int> ();
            _threads = _xtpdft_options.get(key + ".threads").as<int> ();
            _cleanup = _xtpdft_options.get(key + ".cleanup").as<std::string> ();
           
            _write_guess=_xtpdft_options.ifExistsReturnElseReturnDefault<bool>(key + ".read_guess", false);
            
            // check if ECPs are used in xtpdft
            _write_pseudopotentials=false;
            if (_xtpdft_options.exists(key + ".ecp")){
                if (_xtpdft_options.get(key + ".ecp").as<std::string> () !="") {
                    _write_pseudopotentials=true;
                }
            }

        }

        /**
         * Dummy for use of XTPDFT as QMPackage, needs no input file
         */
        bool XTPDFT::WriteInputFile(Orbitals& orbitals) {
            return true;
        }

    
        /**
         * Run calls DFTENGINE
         */
        bool XTPDFT::Run( Orbitals& orbitals ) {
          DFTEngine xtpdft;
          xtpdft.Initialize(_xtpdft_options);
          xtpdft.setLogger(_pLog);
           
          if(_write_charges){
            xtpdft.setExternalcharges(_PolarSegments);
          }
          xtpdft.Prepare( orbitals );
          xtpdft.Evaluate( orbitals );
          _basisset_name = xtpdft.getDFTBasisName();
          orbitals.WriteToCpt(_log_file_name);
          return true;

        }

    void XTPDFT::CleanUp() {
      if (_cleanup.size() != 0) {
        XTP_LOG(logDEBUG, *_pLog) << "Removing " << _cleanup << " files" << flush;
        tools::Tokenizer tok_cleanup(_cleanup, ", ");
        std::vector <std::string> cleanup_info;
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
        bool XTPDFT::ParseOrbitalsFile(Orbitals & orbitals) {
            return true;
        }

        /**
         * Dummy, because information is directly stored in orbitals
         */
        bool XTPDFT::ParseLogFile(Orbitals & orbitals) {
          try{
          orbitals.ReadFromCpt(_log_file_name);
          }catch(std::runtime_error& error){
            XTP_LOG(logDEBUG, *_pLog) << "Reading"<<_log_file_name<<" failed" << flush;
            return false;
          }
            return true;
        }



    }
}
