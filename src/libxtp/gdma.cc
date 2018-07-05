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



#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/algorithm/string.hpp>
#include <votca/ctp/logger.h>
#include <votca/xtp/gdma.h>




namespace votca {
    namespace xtp {
      
      using namespace std;

        // initialize the GDMA object, set parameters
        // for use of external code -> Rank

        void GDMA::Initialize(tools::Property* options) {

            string key = "gdma";
            if (options->exists(key + ".chk")) {
                _chkFile = options->get(key + ".chk").as<string> ();
            } else {
                _chkFile = "system.chk";
            }

            if (options->exists(key + ".executable")) {
                _executable = options->get(key + ".executable").as<string> ();
            } else {
                _executable = "gdma";
            }

            if (options->exists(key + ".density")) {
                _density = options->get(key + ".density").as<string> ();
            } else {
                _density = "SCF"; // defaults to ground state SCF density
            }

            if (options->exists(key + ".limit")) {
                _limit = options->get(key + ".multipoles.limit").as<int> ();
            } else {
                _limit = 2; // default to quadrupoles
            }

            if (_limit > 2) {
                cerr << "Tried to use GDMA with Rank > 2 ";
                throw std::runtime_error("Not supported!");
            }

            // specifying this as a single double seems stupid because it's 
            // element-based...
            // _radius = options->get(key + ".multipoles.radius").as<double> ();

            if (options->exists(key + ".switch")) {
                _switch = options->get(key + ".multipoles.switch").as<double> ();
            } else {
                _switch = 4.0; // corresponds to GDMA default
            }

            if (options->exists(key + ".output")) {
                _outFile = options->get(key + ".output").as<string> ();
            } else {
                _outFile = "gdma.out";
            }

        }

        // write an input file for the external GDMA code by A. Stone

        void GDMA::WriteInputFile() {

            // CTP_LOG(logINFO, *_log) << "Running GDMA " << flush;
            // prepare a GDMA input file
            ofstream _gdma_inputfile;
            string _gdma_inputfile_name_full = _runFolder + "/gdma.in";
            _gdma_inputfile.open(_gdma_inputfile_name_full.c_str());
            _gdma_inputfile << "Title \"Multipole Fit for QMMM\"" << endl;
            _gdma_inputfile << "File system.fchk Density " << _density << endl;
            _gdma_inputfile << "Angstrom" << endl;
            _gdma_inputfile << "Multipoles" << endl;
            _gdma_inputfile << "Limit " << _limit << endl;
            _gdma_inputfile << "Switch " << _switch << endl;
            _gdma_inputfile << "Start" << endl;
            _gdma_inputfile << "Finish" << endl;
            _gdma_inputfile.close();

        }

        void GDMA::RunExternal() {


            // check if the input file exists
            string fullInput = _runFolder + "/gdma.in";
            if (!boost::filesystem::exists(fullInput)) {
                CTP_LOG(ctp::logINFO, *_log) << "GDMA input file has not been found!" << flush;
                throw runtime_error(" GDMA cannot be run! ");
            }

            // check if fchk exists
            string fullFChk = _runFolder + "/system.fchk";
            if (!boost::filesystem::exists(fullFChk)) {
                // try converting Chk to FChk
                string fullChk = _runFolder + "/" + _chkFile;
                if (boost::filesystem::exists(fullChk)) {
                    // use formchk
                    string _command;
                    _command = "cd " + _runFolder + "; formchk " + _chkFile + " system.fchk > /dev/null";
                    if (std::system(_command.c_str())){
                      throw runtime_error("Command "+ _command + "failed");
                    }
                    // check again for fchk
                    if (!boost::filesystem::exists(fullFChk)) {
                        CTP_LOG(ctp::logINFO, *_log) << "Formatted Checkpoint file has not been found and cannot be created!" << flush;
                        throw runtime_error(" GDMA cannot be run! ");

                    }
                } else {
                    CTP_LOG(ctp::logINFO, *_log) << "Formatted Checkpoint file has not been found and cannot be created!" << flush;
                    throw runtime_error(" GDMA cannot be run! ");
                }
            }


            // now we seem ready to go
            string _command;
            _command = "cd " + _runFolder + "; " + _executable + " < gdma.in > " + _outFile;
            if (std::system(_command.c_str())){
              throw runtime_error("Command "+ _command + "failed");
            }

        }

        void GDMA::ParseOutputFile() {

            string _gdma_output_name_full = _runFolder + "/gdma.out";
            std::ifstream _gdma_output(_gdma_output_name_full.c_str());
            std::string _line;
            while (_gdma_output) {

                getline(_gdma_output, _line);
                // if a line has an equality sign, must be energy
                std::string::size_type atom_pos = _line.find("x =");
                if (atom_pos != std::string::npos) {
                    // found an atom, read one more line
                    getline(_gdma_output, _line);
                    // determine rank
                    std::vector<string> results;
                    boost::trim(_line);
                    std::vector<double> Qs; // temp vector for reading in
                    
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "),
                            boost::algorithm::token_compress_on);
                    //cout << results[1] << ":" << results[2] << ":" << results[3] << ":" << results[4] << endl;

                    int rank = boost::lexical_cast<int>(results[3]);
                    if (rank < 0 || rank > 2) {
                        throw runtime_error((boost::format(" Invalid GDMA rank %s!") %rank).str());
                    }

                    // atomic charge
                    if (rank >= 0) {
                        // getting charge
                        getline(_gdma_output, _line);
                        std::vector<string> results;
                        boost::trim(_line);
                        boost::algorithm::split(results, _line, boost::is_any_of("\t "),
                                boost::algorithm::token_compress_on);
                        double Q00 = boost::lexical_cast<double>(results.back());
                        Qs.push_back(Q00);

                        // CTP_LOG(logINFO, *_log) << "New Q00 " << Q00 << flush;

                    }

                    // atomic dipole components
                    if (rank >= 1) {

                        // getting dipoles
                        getline(_gdma_output, _line);
                        std::vector<string> results;
                        boost::trim(_line);
                        boost::algorithm::split(results, _line, boost::is_any_of("\t "),
                                boost::algorithm::token_compress_on);
                        double Q10 = boost::lexical_cast<double>(results[5]);
                        double Q11c = boost::lexical_cast<double>(results[8]);
                        double Q11s = boost::lexical_cast<double>(results.back());
                        Qs.push_back(Q10);
                        Qs.push_back(Q11c);
                        Qs.push_back(Q11s);
                        /*
                        CTP_LOG(logINFO, *_log) << "New Q10  " << Q10 << flush;
                        CTP_LOG(logINFO, *_log) << "New Q11c " << Q11c << flush;
                        CTP_LOG(logINFO, *_log) << "New Q11s " << Q11s << flush;
                         */

                    }
                    
                    // atomic quadrupole components
                    if (rank == 2) {
                        
                        getline(_gdma_output, _line);
                        std::vector<string> results;
                        boost::trim(_line);
                        boost::algorithm::split(results, _line, boost::is_any_of("\t "),
                                boost::algorithm::token_compress_on);
                        double Q20 = boost::lexical_cast<double>(results[5]);
                        double Q21c = boost::lexical_cast<double>(results[8]);
                        double Q21s = boost::lexical_cast<double>(results.back());
                        getline(_gdma_output, _line);
                        boost::trim(_line);
                        boost::algorithm::split(results, _line, boost::is_any_of("\t "),
                                boost::algorithm::token_compress_on);
                        double Q22c = boost::lexical_cast<double>(results[2]);
                        double Q22s = boost::lexical_cast<double>(results[5]);
                        
                        Qs.push_back(Q20);
                        Qs.push_back(Q21c);
                        Qs.push_back(Q21s);
                        Qs.push_back(Q22c);
                        Qs.push_back(Q22s);

                        /* 
                        CTP_LOG(logINFO, *_log) << "New Q20  " << Q20 << flush;
                        CTP_LOG(logINFO, *_log) << "New Q21c " << Q21c << flush;
                        CTP_LOG(logINFO, *_log) << "New Q21s " << Q21s << flush;
                        CTP_LOG(logINFO, *_log) << "New Q22c " << Q22c << flush;
                        CTP_LOG(logINFO, *_log) << "New Q22s " << Q22s << flush;
                        */
                        
                    }

                    _multipoles.push_back( Qs );
                    
                } // atom

                std::string::size_type break_pos = _line.find("Total multipoles referred to origin at");
                if (break_pos != std::string::npos) break;


            } // gdma_output
           // CTP_LOG(logINFO, *_log) << "Done with GDMA" << flush;


        }












    }
}
