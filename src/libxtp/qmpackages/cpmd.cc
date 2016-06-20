/*
 *            Copyright 2009-2016 The VOTCA Development Team
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

#include "cpmd.h"
#include "votca/xtp/segment.h"

#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/constants.h>
#include <stdio.h>
#include <iomanip>
#include <sys/stat.h>
#include <vector>



namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        void Cpmd::Initialize(Property *options) {
            
            //TODO: Yuriy, fill this in
        }

        bool Cpmd::WriteInputFile(std::vector<Segment* > segments, Orbitals* orbitals_guess) {
            
            //TODO: Yuriy, fill this in
            
            return true;
        }

        /**
         * Runs the CPMD job. 
         */
        bool Cpmd::Run() {

            LOG(logDEBUG, *_pLog) << "CPMD: running [" << _executable << " " << _input_file_name << "]" << flush;

            if (std::system(NULL)) {
                
                _command = "cd " + _run_dir + "; mkdir -p $CPMD_SCRDIR; " + _executable + " " + _input_file_name;
                std::system(_command.c_str());
                
                if (CheckLogFile()) {
                    LOG(logDEBUG, *_pLog) << "CPMD: finished job" << flush;
                    return true;
                } else {
                    LOG(logDEBUG, *_pLog) << "CPMD: job failed" << flush;
                }
            } else {
                LOG(logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                return false;
            }

            return true;

        }

        /**
         * Cleans up after the CPMD job
         */
        void Cpmd::CleanUp() {

            //TODO: Yuriy, fill this in
           
        }

        

        bool Cpmd::CheckLogFile() {

            // check if the log file exists
            boost::filesystem::path arg_path;
            char ch;

            std::string _full_name = (arg_path / _run_dir / _log_file_name).c_str();
            ifstream _input_file(_full_name.c_str());

            if (_input_file.fail()) {
                LOG(logERROR, *_pLog) << "CPMD: " << _full_name << " is not found" << flush;
                return false;
            };

            //Use brute force. Search every line for the termination string.
            //It doesn't appear at the very end, like in gaussian
            std::string::size_type self_energy_pos;
            std::string _line;
            do {
                getline(_input_file, _line);
            } while (self_energy_pos!=_line.find("PROGRAM CPMD ENDED AT") || _input_file.eof());
            
            _input_file.close();

            if (self_energy_pos == std::string::npos) {
                LOG(logERROR, *_pLog) << "CPMD: " << _full_name << " is incomplete" << flush;
                return false;
            } else {
                //LOG(logDEBUG,*_pLog) << "CPMD LOG is complete" << flush;
                return true;
            }
        }

        /**
         * Parses the CPMD Log file and stores data in the Orbitals object 
         */
        bool Cpmd::ParseLogFile(Orbitals * _orbitals) {

            //TODO: Yuriy, fill this in
            
            return true;
        }




    }
}
