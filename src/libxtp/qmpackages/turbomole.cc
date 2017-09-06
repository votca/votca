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

#include "turbomole.h"
#include <votca/ctp/segment.h>
#include <votca/tools/globals.h>
#include <votca/xtp/qminterface.h>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <stdio.h>
#include <iomanip>
#include <sys/stat.h>


using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        void Turbomole::Initialize(Property *opt) {

            string key = "package";

            string _name = opt->get(key + ".name").as<string> ();

            if (_name != "turbomole") {
                cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error("Package is not supported.");
            }

            _executable = opt->get(key + ".executable").as<string> ();
            _options = opt->get(key + ".options").as<string> ();
            _scratch_dir = opt->get(key + ".scratch").as<string> ();
            _cleanup = opt->get(key + ".cleanup").as<string> ();

            _input_file_name = "input";
            _log_file_name = _executable + ".log";
            _orb_file_name = "mos";
            _xyz_file_name = "coord";

            _get_charges = false;
            _get_self_energy = false;
            _write_guess = false;

            if (opt->exists(key + ".outputVxc")) {
                _output_Vxc = opt->get(key + "outputVxc").as<bool> ();
            } else _output_Vxc = false;
            if (_output_Vxc) {
                throw std::runtime_error("Sorry " + _name + " does not support Vxc output");
            }

            // check if pseudopotentials are required (NEEDS TO BE ADDED)
            //iop_pos = _options.find("pseudo");
            //if (iop_pos != std::string::npos) {
            //    _write_pseudopotentials = true;
             //   _ecp_name = options->get(key + ".ecp").as<std::string> ();
            //} else {
            _write_pseudopotentials = false;
            //}

            // check if the guess keyword is present, if yes, append the guess later
            std::string::size_type iop_pos1 = _options.find("iter\n1 "); // for 1 + space
            std::string::size_type iop_pos2 = _options.find("iter\n1\n"); // for 1 + new line
            if (iop_pos1 != std::string::npos || iop_pos2 != std::string::npos) _write_guess = true;

        }

        /**
         * Prepares the com file from a vector of segments
         * Appends a guess constructed from monomer orbitals if supplied
         */
        bool Turbomole::WriteInputFile(std::vector< ctp::Segment* > segments, Orbitals* orbitals_guess , std::vector<ctp::PolarSeg*> PolarSegments ) {
            
            /* No background charge writing is implemented. Throw error when
             * getting a non-zero length vector of PolarSegments
             */
            if ( PolarSegments.size() != 0 ){
                throw std::runtime_error("Turbomole cannot be run with multipole background in this version.");
            }
            
            std::vector< ctp::QMAtom* > qmatoms;
            if (_write_charges) {
                qmatoms = orbitals_guess->QMAtoms();
            } else {
                QMMInterface qmmface;
                qmatoms = qmmface.Convert(segments);
                
            }
            
            std::vector< ctp::Atom* > _atoms;
            std::vector< ctp::Atom* > ::iterator ait;
            std::vector< ctp::Segment* >::iterator sit;
            string temp_suffix = "/id";

            double nm2Bohr = tools::conv::nm2bohr;

            CTP_LOG(ctp::logDEBUG, *_pLog) << "TURBOMOLE: Preparing input " << flush;

            std::ofstream _coord_file;

            string _xyz_file_name_full = _run_dir + "/" + _xyz_file_name;

            //cerr << "FILE NAME: " << _com_file_name_full << endl;

            _coord_file.open(_xyz_file_name_full.c_str());

            // header
            _coord_file << "$coord" << endl;

            for (sit = segments.begin(); sit != segments.end(); ++sit) {

                // PBCs are taken care of here

                _atoms = (*sit)-> Atoms();

                temp_suffix = temp_suffix + "_" + boost::lexical_cast<string>((*sit)->getId());

                for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                    if ((*ait)->HasQMPart() == false) {
                        continue;
                    }

                    vec pos = (*ait)->getQMPos();
                    string name = (*ait)->getElement();

                    //fprintf(out, "%2s %4.7f %4.7f %4.7f \n"
                    _coord_file << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX() * nm2Bohr
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY() * nm2Bohr
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ() * nm2Bohr
                            << setw(3) << name.c_str()
                            << endl;
                }
            }

            _coord_file << "$end" << endl;
            _coord_file.close();

            //cerr << _options << flush;

            // "define" does not override the control file but uses it. trash the old version
            std::string file_name = _run_dir + "/control";
            remove(file_name.c_str());

            std::string _command;
            std::string _input_file_name_full = _run_dir + "/" + _input_file_name;

            std::ofstream _input_file;
            _input_file.open(_input_file_name_full.c_str());
            _input_file << "\n" << _options;
            _input_file.close();

            // run "define" which prepares the input
            std::string _input_exe = "define";
            _command = "cd " + _run_dir + "; " + _input_exe + " <  ./" + _input_file_name + " >& " + _input_file_name + ".log";
            int check = std::system(_command.c_str());
            if (check == -1) {
                CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                return false;
            }

            // postprocess the output of define - scratch dir
            //cout <<  "TEMP DIR: " << _scratch_dir + temp_suffix << endl;

            if (_scratch_dir != "") {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "TURBOMOLE: scratch dir " << _scratch_dir + temp_suffix << flush;

                boost::filesystem::create_directories(_scratch_dir + temp_suffix);

                std::ifstream _input_file;
                std::ofstream _temp_input_file;

                std::string _control_file_name_full = _run_dir + "/control";
                std::string _temp_control_file_name_full = _run_dir + "/control.temp";

                _input_file.open(_control_file_name_full.c_str());
                _temp_input_file.open(_temp_control_file_name_full.c_str());

                std::string _line;
                std::string::size_type _pos;
                while (_input_file) {
                    getline(_input_file, _line);

                    _pos = _line.find("$scfdump");
                    if (_pos != std::string::npos) continue;

                    _pos = _line.find("$end");
                    std::string _temp("$TMPDIR " + _scratch_dir + temp_suffix + "\n");
                    if (_pos != std::string::npos) _temp_input_file << _temp;

                    _temp_input_file << _line << endl;
                }

                remove(_control_file_name_full.c_str());
                rename(_temp_control_file_name_full.c_str(), _control_file_name_full.c_str());

            }


            // prepare guess for the orbitals by merging monomer orbitals
            if (_write_guess) {
                if (orbitals_guess == NULL) {
                    throw std::runtime_error("A guess for dimer orbitals has not been prepared.");
                } else {
                    
                    orbitals_guess->QMAtoms() = qmatoms;
                    std::ofstream _orb_file;
                    std::string _orb_file_name_full = _run_dir + "/" + _orb_file_name;
                    _orb_file.open(_orb_file_name_full.c_str());

                    // header
                    _orb_file << "$scfmo    scfconv=1   format(4d20.14)\n#generated by VOTCA\n#\n" << flush;

                    

                    
                    ReorderMOsBack(orbitals_guess);
                    ub::matrix<double>& MOs=orbitals_guess->MOCoefficients();
                    std::vector<int> _sort_index=orbitals_guess->SortEnergies();
                    

                    int level = 1;
                    int ncolumns = 4;

                    for (std::vector< int > ::iterator soi = _sort_index.begin(); soi != _sort_index.end(); ++soi) {

                        double _energy = (orbitals_guess->MOEnergies())[*soi];

                        ub::matrix_row< ub::matrix<double> > mr(MOs, *soi);

                        _orb_file << setw(6) << level << "  a      eigenvalue=" << FortranFormat(_energy) << "   nsaos=" << mr.size() << endl;

                        int column = 1;
                        for (unsigned j = 0; j < mr.size(); ++j) {
                            _orb_file << FortranFormat(mr[j]);
                            if (column == ncolumns) {
                                _orb_file << endl;
                                column = 0;
                            }
                            column++;
                        }

                        level++;
                        if (column != 1) _orb_file << endl;
                    }
                    _orb_file << "$end" << endl;
                }
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << "TURBOMOLE: Finished with input" << flush;
            return true;

        }

        std::string Turbomole::FortranFormat(const double &number) {
            std::stringstream _ssnumber;
            if (number >= 0) _ssnumber << " ";
            _ssnumber << setiosflags(ios::fixed) << setprecision(13) << std::scientific << number;
            std::string _snumber = _ssnumber.str();
            boost::replace_first(_snumber, "e", "D");
            return _snumber;
        }

        /**
         * Runs the TURBOMOLE job.
         */
        bool Turbomole::Run( Orbitals* _orbitals ) {

            CTP_LOG(ctp::logDEBUG, *_pLog) << "TURBOMOLE: Running job [" << _executable << "]" << flush;

            if (std::system(NULL)) {
                // if scratch is provided, run the shell script;
                // otherwise run gaussian directly and rely on global variables
                std::string _command;
                _command = "cd " + _run_dir + "; " + _executable + " >& " + _executable + ".log ";

                //int i = std::system ( _command.c_str() );
                if (std::system(_command.c_str())) {
                    throw runtime_error("Command " + _command + "failed");
                }
                CTP_LOG(ctp::logDEBUG, *_pLog) << "TURBOMOLE: Finished job" << flush;
                return true;
            } else {
                CTP_LOG(ctp::logERROR, *_pLog) << "TURBOMOLE: " << _input_file_name << " failed to start" << flush;
                return false;
            }




        }

        /**
         * Cleans up after the TURBOMOLE job
         */
        void Turbomole::CleanUp() {

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Removing files " << _cleanup << flush;

            // cleaning up the generated files
            if (_cleanup.size() != 0) {
                Tokenizer tok_cleanup(_cleanup, ",");
                std::vector <std::string> _cleanup_info;
                tok_cleanup.ToVector(_cleanup_info);

                std::vector<std::string> ::iterator it;

                for (it = _cleanup_info.begin(); it != _cleanup_info.end(); ++it) {
                    if (*it == "input") {
                        std::string file_name = _run_dir + "/" + _input_file_name;
                        remove(file_name.c_str());
                    }

                    if (*it == "log") {
                        std::string file_name = _run_dir + "/" + _log_file_name;
                        remove(file_name.c_str());
                    }

                    if (*it == "mos") {
                        std::string file_name = _run_dir + "/" + *it;
                        remove(file_name.c_str());
                    }
                }
            }

        }

        /**
         * Reads in the MO coefficients from a TURBOMOLE file
         */
        bool Turbomole::ParseOrbitalsFile(Orbitals* _orbitals) {
            std::map <int, std::vector<double> > _coefficients;
            std::map <int, double> _energies;

            //double _conv_Hrt_eV = 27.21138386;

            std::string _line;
            unsigned _levels = 0;
            unsigned _level = 0;
            unsigned _basis_size = 0;

            path arg_path;
            std::string orbFileName = (arg_path / _run_dir / _orb_file_name).string();
            std::ifstream _input_file(orbFileName.c_str());
            //cout << endl << (_run_dir + "/" + _orb_file_name).c_str();

            if (_input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "File " << _orb_file_name << " with molecular orbitals is not found " << flush;
                return false;
            } else {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Reading MOs from " << _orb_file_name << flush;
            }

            // get the first line with $
            getline(_input_file, _line);
            // skip all comments (lines with #)

            getline(_input_file, _line);
            std::string::size_type hash_pos = _line.find("#");

            while (hash_pos != std::string::npos) {
                getline(_input_file, _line);
                hash_pos = _line.find("#");
            }

            //clog << endl << "Orbital file " << filename << " has "
            //        << nrecords_in_line << " records per line, in D"
            //        << format << " format." << endl;

            while (_input_file) {

                // if a line has an equality sign, must be energy
                std::string::size_type energy_pos = _line.find("=");
                std::string::size_type dollar_pos = _line.find("$");

                if (energy_pos != std::string::npos) {

                    std::vector<std::string> results;
                    boost::trim(_line);

                    boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                            boost::algorithm::token_compress_on);
                    //cout << results[1] << ":" << results[2] << ":" << results[3] << ":" << results[4] << endl;

                    //exit(0);

                    _level = boost::lexical_cast<int>(results.front());
                    boost::replace_first(results[3], "D", "e");
                    _energies[ _level ] = boost::lexical_cast<double>(results[3]);
                    _levels++;

                } else if (dollar_pos == std::string::npos) {

                    while (_line.size() > 1) {
                        std::string _coefficient;
                        _coefficient.assign(_line, 0, 20);
                        boost::trim(_coefficient);
                        boost::replace_first(_coefficient, "D", "e");
                        double coefficient = boost::lexical_cast<double>(_coefficient);
                        _coefficients[ _level ].push_back(coefficient);
                        _line.erase(0, 20);
                    }
                } else if (dollar_pos != std::string::npos) break;

                getline(_input_file, _line);
            }

            // some sanity checks
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Energy levels: " << _levels << flush;

            std::map< int, std::vector<double> >::iterator iter = _coefficients.begin();
            _basis_size = iter->second.size();

            for (iter = _coefficients.begin()++; iter != _coefficients.end(); iter++) {
                if (iter->second.size() != _basis_size) {
                    CTP_LOG(ctp::logERROR, *_pLog) << "Error reading " << _orb_file_name << ". Basis set size change from level to level." << flush;
                    return false;
                }
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Basis set size: " << _basis_size << flush;

            // copying information to the orbitals object
            _orbitals->setBasisSetSize(_basis_size);
            // _orbitals->_has_mo_coefficients = true;
            // _orbitals->_has_mo_energies = true;

            // copying energies to a matrix
            _orbitals->MOEnergies().resize(_levels);
            _level = 1;
            for (size_t i = 0; i < _orbitals->_mo_energies.size(); i++) {
                _orbitals->MOEnergies()[i] = _energies[ _level++ ];
            }

            // copying orbitals to the matrix
            (_orbitals->MOCoefficients()).resize(_levels, _basis_size);
            for (size_t i = 0; i < _orbitals->MOCoefficients().size1(); i++) {
                for (size_t j = 0; j < _orbitals->MOCoefficients().size2(); j++) {
                    _orbitals->MOCoefficients()(i, j) = _coefficients[i + 1][j];
                    //cout << i << " " << j << endl;
                }
            }
            
            


            //cout << _mo_energies << endl;
            //cout << _mo_coefficients << endl;

            // cleanup
            _coefficients.clear();
            _energies.clear();


            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done reading MOs" << flush;
            
            return true;
        }

        /*
         * Checks completeness of the TURBOMOLE log
         */
        bool Turbomole::CheckLogFile() {

            // check if the log file exists
            char ch;
            path arg_path;
            std::string logFileName = (arg_path / _run_dir / _log_file_name).string();
            std::ifstream _input_file(logFileName.c_str());
            //cout << (_run_dir + "/" + _log_file_name).c_str();
            if (_input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "TURBOMOLE: " << _log_file_name << " is not found" << flush;
                return false;
            };

            _input_file.seekg(0, ios_base::end); // go to the EOF

            // get empty lines and end of lines out of the way
            do {
                _input_file.seekg(-2, ios_base::cur);
                _input_file.get(ch);
                //cout << "\nChar: " << ch << endl;
            } while (ch == '\n' || ch == ' ' || ch == '\t' || (int) _input_file.tellg() == -1);

            // get the beginning of the line or the file
            do {
                _input_file.seekg(-2, ios_base::cur);
                _input_file.get(ch);
                //cout << "\nNext Char: " << ch << " TELL G " <<  (int)_input_file.tellg() << endl;
            } while (ch != '\n' && ((int) _input_file.tellg() != -1));

            std::string _line;
            getline(_input_file, _line); // Read the current line
            //cout << "\nResult: " << _line << '\n';     // Display it
            _input_file.close();

            std::string::size_type self_energy_pos = _line.find("ended normally");
            if (self_energy_pos == std::string::npos) {
                CTP_LOG(ctp::logERROR, *_pLog) << "TURBOMOLE: " << _log_file_name << " is incomplete" << flush;
                return false;
            } else {
                //CTP_LOG(logDEBUG,*_pLog) << "Gaussian LOG is complete" << flush;
                return true;
            }
        }

        /**
         * Parses the Turbomole Log file and stores data in the Orbitals object
         * TO DO
         */
        bool Turbomole::ParseLogFile(Orbitals* _orbitals) {

            std::string _line;
            std::vector<std::string> results;
            bool _has_occupied_levels = false;
            bool _has_unoccupied_levels = false;
            bool _has_number_of_electrons = false;
            bool _has_basis_set_size = false;
            bool _has_overlap_matrix = false;
            bool _has_charges = false;
            //bool _has_coordinates = false;
            //bool _has_qm_energy = false;
            bool _has_self_energy = false;

            //int _occupied_levels = 0;
            //int _unoccupied_levels = 0;
            int _number_of_electrons = 0;
            int _basis_set_size = 0;


            CTP_LOG(ctp::logDEBUG, *_pLog) << "TURBOMOLE: Parsing " << _log_file_name << flush;

            // check if LOG file is complete
            if (!CheckLogFile()) return false;
            // save qmpackage name
            //_orbitals->_has_qm_package = true;
            _orbitals->setQMpackage("turbomole");
            _orbitals->setDFTbasis(_basisset_name);


            if (_write_pseudopotentials) {
                _orbitals->setECP(_ecp_name);
            } 

            // Start parsing the file line by line
            path arg_path;
            std::string logFileName = (arg_path / _run_dir / _log_file_name).string();
            std::ifstream _input_file(logFileName.c_str());

            while (_input_file) {

                getline(_input_file, _line);
                boost::trim(_line);


                /*
                 * number of occupied and virtual orbitals
                 * N alpha electrons      M beta electrons
                 */
                std::string::size_type electrons_pos = _line.find("number of occupied orbitals");
                if (electrons_pos != std::string::npos) {
                    //cout << _line << endl;
                    boost::algorithm::split(results, _line, boost::is_any_of(":"), boost::algorithm::token_compress_on);
                    _has_number_of_electrons = true;
                    std::string _oo = results.back();
                    boost::trim(_oo);
                    _number_of_electrons = boost::lexical_cast<int>(_oo);
                    _orbitals->setNumberOfElectrons(_number_of_electrons);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Alpha electrons: " << _number_of_electrons << flush;
                }

                /*
                 * basis set size
                 * N basis functions,  M primitive gaussians,   K cartesian basis functions
                 */
                std::string::size_type basis_pos = _line.find("number of basis functions");
                if (basis_pos != std::string::npos) {
                    //cout << _line << endl;
                    boost::algorithm::split(results, _line, boost::is_any_of(":"), boost::algorithm::token_compress_on);
                    _has_basis_set_size = true;
                    std::string _bf = results.back();
                    boost::trim(_bf);
                    _basis_set_size = boost::lexical_cast<int>(_bf);
                    _orbitals->setBasisSetSize(_basis_set_size);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Basis functions: " << _basis_set_size << flush;
                }

                /*
                 * number of unoccupied orbitals = basis set size - number of occupied orbitals
                 */
                if (_has_basis_set_size && _has_number_of_electrons) {
                    _orbitals->setNumberOfLevels(_number_of_electrons, _basis_set_size - _number_of_electrons);
                }

                /*
                 * overlap matrix
                 * stored after the *** Overlap *** line
                 */
                std::string::size_type overlap_pos = _line.find("OVERLAP");
                if (overlap_pos != std::string::npos) {

                    // prepare the container
                    // _orbitals->_has_overlap = true;
                    (_orbitals->_overlap).resize(_basis_set_size);

                    _has_overlap_matrix = true;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Read the overlap matrix" << flush;

                    // skip the next line with "----"
                    getline(_input_file, _line);
                    getline(_input_file, _line);

                    overlap_pos = _line.find("--");

                    int _i_index = 0;
                    int _j_index = 0;

                    while (overlap_pos == std::string::npos) {


                        boost::trim(_line);

                        std::vector<std::string> _row;
                        boost::algorithm::split(_row, _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
                        //int nfields =  _row.size();
                        //cout << nfields << endl;

                        for (std::vector<std::string>::iterator it = _row.begin(); it < _row.end(); it++) {

                            //cout << "  " << *it << endl;

                            boost::trim(*it);
                            double _coefficient = boost::lexical_cast<double>(*it);
                            //cout << _i_index << ":" << _j_index << ":" << _coefficient << endl;

                            _orbitals->AOOverlap()(_i_index, _j_index) = boost::lexical_cast<double>(_coefficient);

                            _j_index++;
                            if (_j_index > _i_index) {
                                _j_index = 0;
                                _i_index++;
                            }
                        }

                        getline(_input_file, _line);
                        overlap_pos = _line.find("--");
                        
                    }


                } // end of the if "Overlap" found


                /*
                 *  Partial charges from the input file [TO DO]
                 */
                std::string::size_type charge_pos = _line.find("Charges from ESP fit, RMS");

                if (charge_pos != std::string::npos && _get_charges) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting charges" << flush;
                    _has_charges = true;
                    //_orbitals->_has_atoms = true;
                }


                /*
                 * Coordinates of the final configuration [TO DO]
                 */
                std::string::size_type coordinates_pos = _line.find("Test job not archived");

                if (coordinates_pos != std::string::npos) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting the coordinates" << flush;
                    //_has_coordinates = true;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << (boost::format("QM energy[eV]: %4.6f ") % _orbitals->getQMEnergy()).str() << flush;

                    //_orbitals->_has_atoms = true;
                    // _orbitals->_has_qm_energy = true;

                }

                /*
                 * Self-energy of external charges [TO DO]
                 */
                std::string::size_type self_energy_pos = _line.find("Self energy of the charges");

                if (self_energy_pos != std::string::npos) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting the self energy\n";
                    _has_self_energy = true;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Self energy " << _orbitals->getSelfEnergy() << flush;

                    // _orbitals->_has_self_energy = true;
                }

                // check if all information has been accumulated and quit
                if (_has_number_of_electrons &&
                        _has_basis_set_size &&
                        _has_occupied_levels &&
                        _has_unoccupied_levels &&
                        _has_overlap_matrix &&
                        _has_charges &&
                        _has_self_energy
                        ) break;

            } // end of reading the file line-by-line
            ReorderOutput( _orbitals);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "TURBOMOLE: Done parsing" << flush;
            return true;
        }





    }
}
