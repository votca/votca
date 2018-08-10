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

#include "nwchem.h"
#include <votca/ctp/segment.h>
#include <votca/xtp/qminterface.h>
#include <votca/xtp/basisset.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <stdio.h>
#include <iomanip>




namespace votca {
    namespace xtp {
       using namespace std;

        void NWChem::Initialize(tools::Property &options) {

            // NWChem file names
            string fileName = "system";

            _xyz_file_name = fileName + ".xyz";
            _input_file_name = fileName + ".nw";
            _log_file_name = fileName + ".log";
            _shell_file_name = fileName + ".sh";
            _orb_file_name = fileName + ".movecs";

            string key = "package";
            string _name = options.get(key + ".name").as<string> ();

            if (_name != "nwchem") {
                cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error("Wrong options file");
            }

            _executable = options.get(key + ".executable").as<string> ();
            _charge = options.get(key + ".charge").as<int> ();
            _spin = options.get(key + ".spin").as<int> ();
            _options = options.get(key + ".options").as<string> ();
            _memory = options.get(key + ".memory").as<string> ();
            _threads = options.get(key + ".threads").as<int> ();
            _scratch_dir = options.get(key + ".scratch").as<string> ();
            _cleanup = options.get(key + ".cleanup").as<string> ();
            
            _basisset_name = options.get(key + ".basisset").as<std::string> ();
            _write_basis_set = options.get(key + ".writebasisset").as<bool> ();
            _write_pseudopotentials = options.get(key + ".writepseudopotentials").as<bool> ();
            
            if ( _write_pseudopotentials )  _ecp_name = options.get(key + ".ecp").as<std::string> ();

            if (options.exists(key + ".outputVxc")) {
                _output_Vxc = options.get(key + ".outputVxc").as<bool> ();
            } else _output_Vxc = false;
            // check whether options string contains vxc output, the _outputVxc is set to true
            std::string::size_type iop_pos = _options.find(" intermediate tXC matrix");
            if (iop_pos != std::string::npos) {
                if (_output_Vxc) {
                    cout << "=== You do not have to specify outputting Vxc twice. Next time remove "
                            "the print ""intermediate tXC matrix"" part from your options string. Please continue" << endl;
                } else {
                    cout << "=== So you do not want to output Vxc but still put it in the options string? "
                            "I will assume that you want to output Vxc, be more consistent next time. " << endl;
                }
                _output_Vxc = true;
            }
            else if (_output_Vxc == true) {
                _options = _options + "\n\ndft\nprint \"intermediate tXC matrix\"\nvectors input system.movecs\nnoscf\nend\ntask dft";
            }

            // check if the optimize keyword is present, if yes, read updated coords
            iop_pos = _options.find(" optimize");
            if (iop_pos != std::string::npos) {
                _is_optimization = true;
            } else {
                _is_optimization = false;
            }

            // check if the esp keyword is present, if yes, get the charges and save them
            iop_pos = _options.find(" esp");
            if (iop_pos != std::string::npos) {
                _get_charges = true;
            } else {
                _get_charges = false;
            }
  

            // check if the guess should be prepared, if yes, append the guess later
            _write_guess = false;
            iop_pos = _options.find("iterations 1 ");
            if (iop_pos != std::string::npos) _write_guess = true;
            iop_pos = _options.find("iterations 1\n");
            if (iop_pos != std::string::npos) _write_guess = true;
        }

         /* For QM/MM the molecules in the MM environment are represented by
         * their atomic partial charge distributions. Triggered by the option
         * keyword "set bq background" NWChem expects them in x,y,z,q format in the
         * backround.crg file.
         */
        
        void NWChem::WriteChargeOption(){
              std::string::size_type iop_pos = _options.find("set bq background");
              if (iop_pos != std::string::npos) {
                _options = _options + "\n set bq background";
              }
        }
       

        int NWChem::WriteBackgroundCharges(ofstream& nw_file) {

            int numberofcharges=0;
            boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
              for (std::shared_ptr<ctp::PolarSeg> seg:_PolarSegments) {
                for (ctp::APolarSite* site:*seg) {
                    string sitestring=boost::str(fmt % ((site->getPos().getX())*votca::tools::conv::nm2ang) 
                            % (site->getPos().getY()*votca::tools::conv::nm2ang) 
                            % (site->getPos().getZ()*votca::tools::conv::nm2ang) 
                            % site->getQ00());
                    if (site->getQ00() != 0.0){
                      nw_file << sitestring << endl;
                      numberofcharges++;
                    }
                    if (site->getRank() > 0 || _with_polarization ) {
                        std::vector< std::vector<double> > _split_multipoles = SplitMultipoles(site);
                        for (const auto& mpoles:_split_multipoles){
                           string multipole=boost::str( fmt % mpoles[0] % mpoles[1] % mpoles[2] % mpoles[3]);
                            nw_file << multipole << endl;
                            numberofcharges++;

                        }
                    }
                }
            }
            nw_file << endl;
            return numberofcharges;
        }
        

        bool NWChem::WriteGuess(Orbitals& orbitals){
          ofstream orb_file;
          std::string orb_file_name_full = _run_dir + "/" + _orb_file_name;
          // get name of temporary ascii file and open it
          std::vector<std::string> results;
          boost::algorithm::split(results, _orb_file_name, boost::is_any_of("."), boost::algorithm::token_compress_on);
          std::string orb_file_name_ascii = _run_dir + "/" + results.front() + ".mos";
          orb_file.open(orb_file_name_ascii.c_str());
          
          // header
          orb_file << "#generated by VOTCA\nbasisum\ngeomsum\n\nscf\nFri Sep 13 00:00:00 2013\nscf\n1\n\n8\nao basis\n1\n" << flush;
          int size_of_basis = (orbitals.MOEnergies()).size();
          orb_file << size_of_basis << endl;
          orb_file << size_of_basis << endl;
          ReorderMOsBack(orbitals);
          int level = 1;
          int ncolumns = 3;
          // write occupations as double in three columns
          // occupied levels
          int column = 1;
          for (int i = 0; i < orbitals.getNumberOfElectrons(); i++) {
            orb_file << FortranFormat(2.0);
            if (column == ncolumns) {
              orb_file << endl;
              column = 0;
            }
            column++;
          }
          // unoccupied levels
          for (int i = orbitals.getNumberOfElectrons(); i < size_of_basis; i++) {
            orb_file << FortranFormat(0.0);
            if (column == ncolumns) {
              orb_file << endl;
              column = 0;
            }
            column++;
          }
          // extra endl
          if (column != 1) {
            orb_file << endl;
          }
          
          // write all energies in same format
          column = 1;
          for (int i=0;i<orbitals.MOEnergies().size();++i) {
            double energy = (orbitals.MOEnergies())[i];
            orb_file << FortranFormat(energy);
            if (column == ncolumns) {
              orb_file << endl;
              column = 0;
            }
            column++;
          }
          if (column != 1) orb_file << endl;
          
          // write coefficients in same format
          for (int i=0;i<orbitals.MOCoefficients().cols();++i) {
            Eigen::VectorXd mr=orbitals.MOCoefficients().col(i);
            column = 1;
            for (unsigned j = 0; j < mr.size(); ++j) {
              orb_file << FortranFormat(mr[j]);
              if (column == ncolumns) {
                orb_file << endl;
                column = 0;
              }
              column++;
            }
            level++;
            if (column != 1) orb_file << endl;
          }
          orb_file << " 0.0000   0.0000" << endl;
          orb_file.close();
          // now convert this ascii file to binary
          std::string command;
          command = "cd " + _run_dir + "; asc2mov 5000 system.mos system.movecs > convert.log";
          int i = std::system(command.c_str());
          if (i == 0) {
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Converted MO file from ascii to binary" << flush;
          } else {
            CTP_LOG(ctp::logERROR, *_pLog) << "Conversion of binary MO file to binary failed. " << flush;
            return false;
          }
          return true;
        }

        /**
         * Prepares the *.nw file from a vector of segments
         * Appends a guess constructed from monomer orbitals if supplied
         */
        bool NWChem::WriteInputFile(Orbitals& orbitals){

           
            std::string temp_suffix = "/id";
            std::string scratch_dir_backup = _scratch_dir;
            std::ofstream nw_file;
            std::ofstream crg_file;

            std::string nw_file_name_full = _run_dir + "/" + _input_file_name;
            std::string crg_file_name_full = _run_dir + "/background.crg";

            nw_file.open(nw_file_name_full.c_str());
            // header
            nw_file << "geometry noautoz noautosym" << endl;

            std::vector< QMAtom* > qmatoms = orbitals.QMAtoms();
          
            for (const QMAtom* atom:qmatoms) {
                tools::vec pos=atom->getPos()*tools::conv::bohr2ang;
                    nw_file << setw(3) << atom->getType().c_str()
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX()
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY()
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ()
                            << endl;      
            }
            nw_file << "end\n";
            if (_write_charges) {
                // part for the MM charge coordinates
                crg_file.open(crg_file_name_full.c_str());
                int numberofcharges=WriteBackgroundCharges(crg_file);
                crg_file << endl;
                crg_file.close();
                nw_file<<endl;
                nw_file<<"set bq:max_nbq "<<numberofcharges<<endl;
                nw_file<<"bq background"<<endl;
                nw_file<<"load background.crg format 1 2 3 4"<<endl;
                nw_file<<"end\n"<<endl;
            }
            
            if(_write_basis_set){
              WriteBasisset(nw_file,qmatoms);
            }
            
            if(_write_pseudopotentials){
              WriteECP(nw_file,qmatoms);
            }

            // write charge of the molecule
            nw_file << "\ncharge " << _charge << "\n";
           
            // writing scratch_dir info
            if (_scratch_dir != "") {
                std::string _temp("scratch_dir " + _scratch_dir + temp_suffix + "\n");
                nw_file << _temp;
            }
            if(_charge!=0.0){
              std::string dft="dft";
              if(_options.find(dft) != std::string::npos){
                int dftpos=_options.find(dft);
                dftpos+=dft.size();
                std::string openshell="\nodft\n" +(boost::format("mult %1%\n") % _spin).str();
                _options.insert(dftpos,openshell,0,openshell.size());
              }else{
                throw runtime_error("NWCHEM: dft input data missing");     
              }
            }
            nw_file << _options << "\n";
            if (_write_guess) {         
                 bool worked=WriteGuess(orbitals);
                 if(!worked){
                    return false;
                 }   
            }

            nw_file << endl;
            nw_file.close();
            // and now generate a shell script to run both jobs, if neccessary
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

            _scratch_dir = scratch_dir_backup + temp_suffix;

            WriteShellScript();
            _scratch_dir = scratch_dir_backup;

            return true;
        }
        

        

        bool NWChem::WriteShellScript() {
            ofstream shell_file;
            std::string shell_file_name_full = _run_dir + "/" + _shell_file_name;
            shell_file.open(shell_file_name_full.c_str());
            shell_file << "#!/bin/bash" << endl;
            shell_file << "mkdir -p " << _scratch_dir << endl;
            if (_threads == 1) {
                shell_file << _executable << " " << _input_file_name << " > " << _log_file_name << " 2> run.error" << endl;
            } else {
                shell_file << "mpirun -np " << boost::lexical_cast<std::string>(_threads) << " " << _executable << " " << _input_file_name << " > " << _log_file_name << " 2> run.error" << endl;
            }
            shell_file.close();
            return true;
        }

        /**
         * Runs the NWChem job.
         */
        bool NWChem::Run( Orbitals& orbitals ) {

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Running NWChem job" << flush;

            if (std::system(NULL)) {

                // NWChem overrides input information, if *.db and *.movecs files are present
                // better trash the old version
                std::string file_name = _run_dir + "/system.db";
                remove(file_name.c_str());
                file_name = _run_dir + "/" + _log_file_name;
                remove(file_name.c_str());
               
                std::string command = "cd " + _run_dir + "; sh " + _shell_file_name;
    
                int check = std::system(command.c_str());
                if (check == -1) {
                    CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                    return false;
                }

                if (CheckLogFile()) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Finished NWChem job" << flush;
                    return true;
                } else {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "NWChem job failed" << flush;
                    return false;
                }
            } else {
                CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                return false;
            }

            return true;
        }

        /**
         * Cleans up after the NWChem job
         */
        void NWChem::CleanUp() {

            // cleaning up the generated files
            if (_cleanup.size() != 0) {
                tools::Tokenizer tok_cleanup(_cleanup, ",");
                std::vector <std::string> cleanup_info;
                tok_cleanup.ToVector(cleanup_info);

                std::vector<std::string> ::iterator it;

                for (const std::string& substring:cleanup_info) {
                    if (substring== "nw") {
                        std::string file_name = _run_dir + "/" + _input_file_name;
                        remove(file_name.c_str());
                    }

                    if (substring == "db") {
                        std::string file_name = _run_dir + "/system.db";
                        remove(file_name.c_str());
                    }

                    if (substring == "log") {
                        std::string file_name = _run_dir + "/" + _log_file_name;
                        remove(file_name.c_str());
                    }

                    if (substring == "movecs") {
                        std::string file_name = _run_dir + "/" + _orb_file_name;
                        remove(file_name.c_str());
                    }

                    if (substring == "gridpts") {
                        std::string file_name = _run_dir + "/system.gridpts.*";
                        remove(file_name.c_str());
                    }
                }
            }

        }

        /**
         * Reads in the MO coefficients from a NWChem movecs file
         */
        bool NWChem::ParseOrbitalsFile(Orbitals& orbitals) {
            std::map <int, std::vector<double> > coefficients;
            std::map <int, double> energies;
            std::map <int, double> occupancy;

            std::string line;
            unsigned levels = 0;
            //unsigned _level;
            unsigned basis_size = 0;
            int number_of_electrons = 0;

            /* maybe we DO need to convert from fortran binary to ASCII first to avoid
             compiler-dependent issues */
            std::string orb_file_name_bin = _run_dir + "/" + _orb_file_name;
            std::string command;
            command = "cd " + _run_dir + "; mov2asc 10000 system.movecs system.mos > convert.log";
            int i = std::system(command.c_str());
            if (i == 0) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Converted MO file from binary to ascii" << flush;
            } else {
                CTP_LOG(ctp::logERROR, *_pLog) << "Conversion of binary MO file to ascii failed. " << flush;
                return false;
            }

            // opening the ascii MO file
            std::string orb_file_name_full = _run_dir + "/" + "system.mos";
            std::ifstream input_file(orb_file_name_full.c_str());

            if (input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "File " << orb_file_name_full << " with molecular orbitals is not found " << flush;
                return false;
            } else {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Reading MOs from " << orb_file_name_full << flush;
            }

            // the first 12 lines are garbage info
            for (i = 1; i < 13; i++) {
                getline(input_file, line);
            }
            // next line has basis set size
            input_file >> basis_size;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Basis set size: " << basis_size << flush;


            // next line has number of stored MOs
            input_file >> levels;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Energy levels: " << levels << flush;

            /* next lines contain information about occupation of the MOs
             *  - each line has 3 numbers
             *  - from _occ we can determine the number of electrons/2 */
            int n_lines = ((levels - 1) / 3);
            int n_rest = levels - 3 * n_lines;
            // read in the data
            int imo = 0;
            for (i = 1; i <= n_lines; i++) {
                for (int j = 0; j < 3; j++) {
                    input_file >> occupancy[ imo ];
                    if (occupancy[ imo ] == 2.0) {
                        number_of_electrons++;
                    }
                    imo++;
                }
            }
            if (n_rest != 0) {
                for (i = 0; i < n_rest; i++) {
                    input_file >> occupancy[ imo ];
                    imo++;
                }
            }
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Alpha electrons: " << number_of_electrons << flush;

            int occupied_levels = number_of_electrons;
            int unoccupied_levels = levels - occupied_levels;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Occupied levels: " << occupied_levels << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Unoccupied levels: " << unoccupied_levels << flush;

            // reset index and read MO energies the same way
            imo = 0;
            for (i = 1; i <= n_lines; i++) {
                for (int j = 0; j < 3; j++) {
                    input_file >> energies[ imo ];
                    imo++;
                }
            }
            if (n_rest != 0) {
                for (i = 0; i < n_rest; i++) {
                    input_file >> energies[ imo ];
                    imo++;
                }
            }


            // Now, the same for the coefficients
            double coef;
            for (unsigned imo = 0; imo < levels; imo++) {
                for (i = 1; i <= n_lines; i++) {
                    for (int j = 0; j < 3; j++) {
                        input_file >> coef;
                        coefficients[ imo ].push_back(coef);
                    }
                }
                if (n_rest != 0) {
                    for (i = 0; i < n_rest; i++) {
                        input_file >> coef;
                        coefficients[ imo ].push_back(coef);
                    }
                }
            }




            // copying information to the orbitals object
            orbitals.setBasisSetSize(basis_size);
            orbitals.setNumberOfElectrons(number_of_electrons);
            orbitals.setNumberOfLevels(occupied_levels, unoccupied_levels);
            // copying energies to a matrix
            orbitals.MOEnergies().resize(levels);
            //_level = 1;
            for (int i = 0; i < orbitals.MOEnergies().size(); i++) {
                orbitals.MOEnergies()[i] = energies[ i ];
            }


            // copying orbitals to the matrix
            (orbitals.MOCoefficients()).resize(levels, basis_size);
            for (int i = 0; i < orbitals.MOCoefficients().rows(); i++) {
                for (int j = 0; j < orbitals.MOCoefficients().cols(); j++) {
                    orbitals.MOCoefficients()(j, i) = coefficients[i][j];
                }
            }
            
            // when all is done, we can trash the ascii file
            std::string file_name = _run_dir + "/system.mos";
            remove(file_name.c_str());


            ReorderOutput(orbitals);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done reading MOs" << flush;

            return true;
        }

        bool NWChem::CheckLogFile() {

            // check if the log file exists
            
            ifstream input_file((_run_dir + "/" + _log_file_name).c_str());

            if (input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "NWChem LOG is not found" << flush;
                return false;
            };

            if (input_file.peek() == std::ifstream::traits_type::eof()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "NWChem run failed. Check OpenMPI version!" << flush;
                return false;
            };



            /* Checking the log file is a pain in the *** since NWChem throws an error
             * for our 'iterations 1'  runs (even though it outputs the required data
             * correctly. The only way that works for both scf and noscf runs is to
             * check for "Total DFT energy" near the end of the log file.
             */

            input_file.seekg(0, ios_base::end); // go to the EOF
            char ch;
            std::string::size_type total_energy_pos = std::string::npos;
            std::string::size_type diis_pos = std::string::npos;
            do {
                // get the beginning of the line
                do {
                    input_file.seekg(-2, ios_base::cur);
                    input_file.get(ch);
                    //cout << "\nNext Char: " << ch << " TELL G " <<  (int)_input_file.tellg() << endl;
                } while (ch != '\n' || (int) input_file.tellg() == -1);

                std::string line;
                getline(input_file, line); // Read the current line
                total_energy_pos = line.find("Total DFT energy");
                diis_pos = line.find("diis");
                // whatever is found first, determines the completeness of the file
                if (total_energy_pos != std::string::npos) {
                    return true;
                } else if (diis_pos != std::string::npos) {
                    CTP_LOG(ctp::logERROR, *_pLog) << "NWChem LOG is incomplete" << flush;
                    return false;
                } else {
                    // go to previous line
                    //_input_file.get(ch);
                    do {
                        input_file.seekg(-2, ios_base::cur);
                        input_file.get(ch);
                        //cout << "\nChar: " << ch << endl;
                    } while (ch != '\n' || (int) input_file.tellg() == -1);
                }
            } while (total_energy_pos == std::string::npos && diis_pos == std::string::npos);


            input_file.close();
            return true;
        }

        

        /**
         * Parses the Gaussian Log file and stores data in the Orbitals object
         */
        bool NWChem::ParseLogFile(Orbitals& orbitals) {

            double conv_Hrt_eV = tools::conv::hrt2ev;

            std::string line;
            std::vector<std::string> results;

            bool has_overlap_matrix = false;
            bool has_charges = false;
            bool has_qm_energy = false;
            bool has_self_energy = false;
            bool has_basis_set_size = false;

            bool found_optimization = false;
            int basis_set_size = 0;

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Parsing " << _log_file_name << flush;
            std::string log_file_name_full = _run_dir + "/" + _log_file_name;
            // check if LOG file is complete
            if (!CheckLogFile()) return false;

            // save qmpackage name
            orbitals.setQMpackage("nwchem");
            orbitals.setDFTbasis(_basisset_name);



            if (_write_pseudopotentials) {
                orbitals.setECP(_ecp_name);
            } 
            // set _found_optimization to true if this is a run without optimization
            if (!_is_optimization) {
                found_optimization = true;
            }

            // Start parsing the file line by line
            ifstream input_file(log_file_name_full.c_str());
            while (input_file) {

                getline(input_file, line);
                boost::trim(line);

                /*
                 * basis set size (is only required for overlap matrix reading, rest is
                 * in orbitals file and could be skipped
                 */
                std::string::size_type basis_pos = line.find("number of functions");
                if (basis_pos != std::string::npos) {
                    //cout << _line << endl;
                    boost::algorithm::split(results, line, boost::is_any_of(":"), boost::algorithm::token_compress_on);
                    has_basis_set_size = true;
                    std::string bf = results.back();
                    boost::trim(bf);
                    basis_set_size = boost::lexical_cast<int>(bf);
                    orbitals.setBasisSetSize(basis_set_size);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Basis functions: " << basis_set_size << flush;
                }

                /*
                 * Total DFT energy
                 */
                std::string::size_type energy_pos = line.find("Total DFT energy");
                if (energy_pos != std::string::npos) {
                    boost::algorithm::split(results, line, boost::is_any_of("="), boost::algorithm::token_compress_on);
                    std::string energy = results.back();
                    boost::trim(energy);
                    orbitals.setQMEnergy(conv_Hrt_eV * boost::lexical_cast<double>(energy));
                    CTP_LOG(ctp::logDEBUG, *_pLog) << (boost::format("QM energy[eV]: %4.6f ") % orbitals.getQMEnergy()).str() << flush;
                    has_qm_energy = true;
                    // _orbitals._has_qm_energy = true;

                }

                /*
                 *  Partial charges from the input file
                 */
                std::string::size_type charge_pos = line.find("ESP");
                if (charge_pos != std::string::npos && _get_charges) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting charges" << flush;
                    has_charges = true;
                    // two empty lines
                    getline(input_file, line);
                    getline(input_file, line);

                    // now starts the data in format
                    // _id type x y z q

                    std::vector<std::string> row=GetLineAndSplit(input_file, "\t ");
                    int nfields = row.size();

                    while (nfields == 6) {
                        int atom_id = boost::lexical_cast< int >(row.at(0));
                        std::string atom_type = row.at(1);
                        double atom_charge = boost::lexical_cast< double >(row.at(5));
                        row=GetLineAndSplit(input_file, "\t ");
                        nfields = row.size();
                        QMAtom* pAtom;
                        if (orbitals.hasQMAtoms() == false) {
                            pAtom =orbitals.AddAtom(atom_id - 1,atom_type, 0, 0, 0);
                        } else {
                            pAtom = orbitals.QMAtoms().at(atom_id - 1);
                        }
                        pAtom->setPartialcharge(atom_charge);
                        }
                }


                /*
                 * Coordinates of the final configuration
                 * depending on whether it is an optimization or not
                 */


                if (_is_optimization) {
                    std::string::size_type optimize_pos = line.find("Optimization converged");
                    if (optimize_pos != std::string::npos) {
                        found_optimization = true;
                    }
                }
                

                std::string::size_type coordinates_pos = line.find("Output coordinates");

                if (found_optimization && coordinates_pos != std::string::npos) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting the coordinates" << flush;

                    //_has_coordinates = true;
                    bool has_QMAtoms = orbitals.hasQMAtoms();

                    // three garbage lines
                    getline(input_file, line);
                    getline(input_file, line);
                    getline(input_file, line);
                    // now starts the data in format
                    // _id type Qnuc x y z
                     std::vector<std::string> row=GetLineAndSplit(input_file, "\t ");
                    int nfields = row.size();
                    
                    while (nfields == 6) {
                        int atom_id = boost::lexical_cast< int >(row.at(0))-1;
                        std::string _atom_type = row.at(1);
                        double x = boost::lexical_cast<double>(row.at(3));
                        double y = boost::lexical_cast<double>(row.at(4));
                        double z = boost::lexical_cast<double>(row.at(5));
                        tools::vec pos=tools::vec(x,y,z);
                        pos*=tools::conv::ang2bohr;
                        if (has_QMAtoms == false) {
                            orbitals.AddAtom(atom_id,_atom_type, pos);
                        } else{
                            QMAtom* pAtom = orbitals.QMAtoms().at(atom_id);
                            pAtom->setPos(pos); 
                        }
                        atom_id++;
                        row=GetLineAndSplit(input_file, "\t ");
                        nfields = row.size();
                    }
                }        
                
                /*
                 * Vxc matrix
                 * stored after the global array: g vxc
                 */
                if(_output_Vxc){
                std::string::size_type vxc_pos = line.find("global array: g vxc");
                if (vxc_pos != std::string::npos) {


                    // prepare the container
                    Eigen::MatrixXd vxc = orbitals.AOVxc();
                    vxc.resize(basis_set_size,basis_set_size);


                    //_has_vxc_matrix = true;
                    std::vector<int> j_indeces;

                    int n_blocks = 1 + ((basis_set_size - 1) / 6);
                    //cout << _n_blocks;

                    for (int block = 0; block < n_blocks; block++) {
                        // first line is garbage
                        getline(input_file, line);
                        // second line gives the j index in the matrix
                        getline(input_file, line);
                        boost::tokenizer<> tok(line);

                        /// COMPILATION IS BROKEN DUE TO A BUG IN BOOST 1.53
                        std::transform(tok.begin(), tok.end(), std::back_inserter(j_indeces), &boost::lexical_cast<int, std::string>);

                        // third line is garbage again
                        getline(input_file, line);

                        // read the block of max _basis_size lines + the following header
                        for (int i = 0; i < basis_set_size; i++) {
                            std::vector<std::string> row=GetLineAndSplit(input_file, "\t ");
                            int i_index = boost::lexical_cast<int>(row.front());
                            row.erase(row.begin());

                            std::vector<int>::iterator j_iter = j_indeces.begin();

                            for (std::string& coefficient: row) {

                                int j_index = *j_iter;
                                vxc(i_index - 1, j_index - 1) = boost::lexical_cast<double>(coefficient);
                                vxc(j_index - 1, i_index - 1) = boost::lexical_cast<double>(coefficient);
                                j_iter++;

                            }


                        }

                        // clear the index for the next block
                        j_indeces.clear();
                    } // end of the blocks
                    
                    
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Read the Vxc matrix" << flush;

                }
                }

                /* Check for ScaHFX = factor of HF exchange included in functional */
                std::string::size_type HFX_pos = line.find("Hartree-Fock (Exact) Exchange");
                if (HFX_pos != std::string::npos) {
                    boost::algorithm::split(results, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    double ScaHFX = boost::lexical_cast<double>(results.back());
                    orbitals.setScaHFX(ScaHFX);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "DFT with " << ScaHFX << " of HF exchange!" << flush;
                }




                /*
                 * overlap matrix
                 * stored after the global array: Temp Over line
                 */
                std::string::size_type overlap_pos = line.find("global array: Temp Over");
                if (overlap_pos != std::string::npos) {

                    // prepare the container
                    (orbitals.AOOverlap()).resize(basis_set_size,basis_set_size);

                    has_overlap_matrix = true;
                    std::vector<int> j_indeces;

                    int n_blocks = 1 + ((basis_set_size - 1) / 6);

                    for (int block = 0; block < n_blocks; block++) {
                        // first line is garbage
                        getline(input_file, line);
                        // second line gives the j index in the matrix
                        getline(input_file, line);
                        boost::tokenizer<> tok(line);

                        /// COMPILATION IS BROKEN DUE TO A BUG IN BOOST 1.53
                        std::transform(tok.begin(), tok.end(), std::back_inserter(j_indeces), &boost::lexical_cast<int, std::string>);

                        // third line is garbage again
                        getline(input_file, line);

                        // read the block of max _basis_size lines + the following header
                        for (int i = 0; i < basis_set_size; i++) {
                            std::vector<std::string> row=GetLineAndSplit(input_file, "\t ");


                            int i_index = boost::lexical_cast<int>(row.front());
                            row.erase(row.begin());

                            std::vector<int>::iterator j_iter = j_indeces.begin();

                           for (std::string& coefficient: row) {
                                int j_index = *j_iter;
                                orbitals.AOOverlap()(i_index - 1, j_index - 1) = boost::lexical_cast<double>(coefficient);
                                orbitals.AOOverlap()(j_index - 1, i_index - 1) = boost::lexical_cast<double>(coefficient);
                                j_iter++;
                            }


                        }

                        // clear the index for the next block
                        j_indeces.clear();
                    } // end of the blocks
                    
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Read the overlap matrix" << flush;
                } // end of the if "Overlap" found

                /*
                 * TODO Self-energy of external charges
                 */
                std::string::size_type self_energy_pos = line.find("Self energy of the charges");

                if (self_energy_pos != std::string::npos) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting the self energy\n";
                    std::vector<std::string> block;
                    std::vector<std::string> energy;
                    boost::algorithm::split(block, line, boost::is_any_of("="), boost::algorithm::token_compress_on);
                    boost::algorithm::split(energy, block[1], boost::is_any_of("\t "), boost::algorithm::token_compress_on);

                    orbitals.setSelfEnergy(conv_Hrt_eV * boost::lexical_cast<double> (energy[1]));

                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Self energy " << orbitals.getSelfEnergy() << flush;

                }

                // check if all information has been accumulated and quit
                if (has_basis_set_size &&
                        has_overlap_matrix &&
                        has_charges &&
                        has_qm_energy &&
                        has_self_energy
                        ) break;

            } // end of reading the file line-by-line
            
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
            return true;
        }
        
        void NWChem::WriteBasisset(ofstream& nw_file, std::vector<QMAtom*>& qmatoms) {

      std::vector<std::string> UniqueElements = FindUniqueElements(qmatoms);
      BasisSet bs;
      bs.LoadBasisSet(_basisset_name);
      CTP_LOG(ctp::logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;
      nw_file << "basis spherical" << endl;
      for (const std::string& element_name : UniqueElements) {
        const Element& element = bs.getElement(element_name);
        for (const Shell& shell : element) {
          //nwchem can only use S,P,SP,D,F,G shells so we split them up if not SP
          if (!shell.isCombined()) {
            // shell type, number primitives, scale factor
            nw_file << element_name << " " << boost::algorithm::to_lower_copy(shell.getType()) << endl;
            for (const GaussianPrimitive& gaussian : shell) {
              for (unsigned _icontr = 0; _icontr < gaussian._contraction.size(); _icontr++) {
                if (gaussian._contraction[_icontr] != 0.0) {
                  nw_file << FortranFormat(gaussian._decay) << " " << FortranFormat(gaussian._contraction[_icontr]) << endl;
                }
              }
            }
          } else {
            string type = shell.getType();
            for (unsigned i = 0; i < type.size(); ++i) {
              string subtype = string(type, i, 1);
              nw_file << element_name << " " << boost::algorithm::to_lower_copy(subtype) << endl;

              for (const GaussianPrimitive& gaussian : shell) {
                nw_file << FortranFormat(gaussian._decay) << " " << FortranFormat(gaussian._contraction[FindLmax(subtype)]) << endl;
              }
            }
          }
        }
      }
      nw_file << "end\n";
      nw_file << endl;

      return;
    }
        
    void NWChem::WriteECP(ofstream& nw_file, std::vector<QMAtom*>& qmatoms) {

      std::vector<std::string> UniqueElements = FindUniqueElements(qmatoms);

      BasisSet ecp;
      ecp.LoadPseudopotentialSet(_ecp_name);

      CTP_LOG(ctp::logDEBUG, *_pLog) << "Loaded Pseudopotentials " << _ecp_name << flush;

      for (const std::string& element_name : UniqueElements) {
        try {
          const Element& element = ecp.getElement(element_name);
        } catch (std::runtime_error& error) {
          CTP_LOG(ctp::logDEBUG, *_pLog) << "No pseudopotential for " << element_name << " available" << flush;
          continue;
        }
        const Element& element = ecp.getElement(element_name);
        // element name, [possibly indeces of centers], zero to indicate the end
        nw_file << element_name << " nelec " << element.getNcore() << endl;
        for (const Shell& shell : element) {
          string shelltype = shell.getType();
          if (shell.getLmax() == element.getLmax()) {
            shelltype = "ul";
          }
          nw_file << element_name << " " << shelltype << endl;
          for (const GaussianPrimitive& gaussian : shell) {
            nw_file << "    " << gaussian._power << " " << FortranFormat(gaussian._decay) << " " << FortranFormat(gaussian._contraction[0]) << endl;
          }
        }
      }
      nw_file << "end\n";
      nw_file << endl;
      return;
    }

             

        std::string NWChem::FortranFormat(const double &number) {
            std::stringstream ssnumber;
            if (number >= 0) {
                ssnumber << "    ";
            } else {
                ssnumber << "   ";
            }
            ssnumber << setiosflags(ios::fixed) << setprecision(15) << std::scientific << number;
            std::string snumber = ssnumber.str();
            return snumber;
        }




    }
}
