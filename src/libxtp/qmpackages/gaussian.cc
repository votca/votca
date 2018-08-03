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

#include "gaussian.h"
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

        void Gaussian::Initialize(tools::Property &options) {

            // GAUSSIAN file names
            std::string fileName = "system";

            _xyz_file_name = fileName + ".xyz";
            _input_file_name = fileName + ".com";
            _log_file_name = fileName + ".log";
            _shell_file_name = fileName + ".sh";
            _orb_file_name = "fort.7";
            _input_vxc_file_name = fileName + "-2.com";


            std::string key = "package";
            std::string _name = options.get(key + ".name").as<std::string> ();

            if (_name != "gaussian") {
                cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error("Wrong options file");
            }

            _executable = options.get(key + ".executable").as<std::string> ();
            _charge = options.get(key + ".charge").as<int> ();
            _spin = options.get(key + ".spin").as<int> ();
            _options = options.get(key + ".options").as<std::string> ();
            _memory = options.get(key + ".memory").as<std::string> ();
            _threads = options.get(key + ".threads").as<int> ();
            _chk_file_name = options.get(key + ".checkpoint").as<std::string> ();
            _scratch_dir = options.get(key + ".scratch").as<std::string> ();
            _cleanup = options.get(key + ".cleanup").as<std::string> ();


            if (options.exists(key + ".vdWRadii")) {
                _vdWfooter = options.get(key + ".vdWRadii").as<std::string> ();
            } else _vdWfooter = "";


            if (options.exists(key + ".outputVxc")) {
                _output_Vxc = options.get(key + ".outputVxc").as<bool> ();
            } else _output_Vxc = false;
            

            /* G09 by default deletes functions from the basisset according to some
             * criterion based on, a.o., the contraction coefficients. This can lead
             * to inconsistencies when MOs are later used in VOTCA's GWBSE modules
             * (and other post-processing routines). G09's default can be modified
             * by the keywork int=nobasistransform. This will add this keyword
             * automatically to the _options string for runs with G09.
             */
            if ( _executable == "g09" ){
                std::string::size_type basistransform_pos = (boost::algorithm::to_lower_copy(_options)).find("nobasistransform");
                if ( basistransform_pos == std::string::npos ){
                    _options = _options + " int=nobasistransform ";
                }
            }
            
            


            // check if the guess keyword is present, if yes, append the guess later
            std::string::size_type iop_pos = _options.find("cards");
            if (iop_pos != std::string::npos) {
                _write_guess = true;
            } else {
                _write_guess = false;
            }

            // check if the pop keyword is present, if yes, get the charges and save them
            iop_pos = _options.find("pop");
            if (iop_pos != std::string::npos) {
                _get_charges = true;
            } else {
                _get_charges = false;
            }

            // check if the charge keyword is present, if yes, get the self energy and save it
            

            // check if the basis set is available ("/gen")
            iop_pos = _options.find("gen");
            if (iop_pos != std::string::npos) {
                _write_basis_set = true;
                _basisset_name = options.get(key + ".basisset").as<std::string> ();
            } else {
                _write_basis_set = false;
            }

            // check if pseudopotentials are required ("pseudo")
            iop_pos = _options.find("pseudo");
            if (iop_pos != std::string::npos) {
                _write_pseudopotentials = true;
                _ecp_name = options.get(key + ".ecp").as<std::string> ();
            } else {
                _write_pseudopotentials = false;
            }

        }

        void Gaussian::WriteChargeOption() {
          std::string::size_type iop_pos = _options.find("charge");
          if (iop_pos == std::string::npos) {
            std::string::size_type pos = _options.find('\n');
            if (pos != std::string::npos) {
              _options.insert(pos, " charge");
            } else {
              _options = _options + " charge";
            }
          }

        }

        /* Custom basis sets are written on a per-element basis to
         * 'elementname'.gbs files, which are then included in the
         * Gaussian input file using @'elementname'.gbs
         */
        void Gaussian::WriteBasisset(std::ofstream& com_file, std::vector<QMAtom*>& qmatoms) {


          std::vector<std::string> UniqueElements= FindUniqueElements(qmatoms);
            BasisSet bs;
            bs.LoadBasisSet(_basisset_name);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;

            for (const std::string& element_name:UniqueElements) {
               
                        const Element& element = bs.getElement(element_name);
                        /* Write each basis set to a element_name.gbs file
                         * and include the gbs file in the com-file via Gaussian's @ function
                         */
                        std::ofstream el_file;
                        std::string el_file_name = _run_dir + "/" + element_name + ".gbs";

                        el_file.open(el_file_name.c_str());
                        // element name, [possibly indeces of centers], zero to indicate the end    
                        com_file << "@" << element_name << ".gbs" << endl;
                        el_file << element_name << " 0" << endl;
                        for (const Shell& shell:element) {
                            //gaussian can only use S,P,SP,D,F,G shells so we split them up if not SP
                            if (shell.getType() == "SP" || !shell.isCombined()) {
                                // shell type, number primitives, scale factor
                                el_file << shell.getType() << " " << shell.getSize() << " " << FortranFormat(shell.getScale()) << endl;
                                for (const GaussianPrimitive& gaussian:shell) {
                                    el_file << FortranFormat(gaussian._decay);
                                    for (const double& contraction:gaussian._contraction) {
                                        if (contraction!= 0.0) {
                                            el_file << " " << FortranFormat(contraction);
                                        }
                                    }
                                    el_file << endl;
                                }                              
                            } else {
                                string type = shell.getType();
                                for (unsigned i = 0; i < type.size(); ++i) {
                                    string subtype = string(type, i, 1);
                                    el_file << subtype << " " << shell.getSize() << " " << FortranFormat(shell.getScale()) << endl;

                                    for (const GaussianPrimitive& gaussian:shell) {
                                        el_file << FortranFormat(gaussian._decay);
                                        el_file << " " << FortranFormat(gaussian._contraction[FindLmax(subtype)]);
                                    }
                                    el_file << endl;  
                                }
                            }
                        }

                        el_file << "****\n";
                        el_file.close();

                    }
                
            com_file << endl;
            return;
        }

        /* If custom ECPs are used, they need to be specified in the input file
         * in a section following the basis set includes.
         */
        void Gaussian::WriteECP(std::ofstream& com_file, std::vector<QMAtom*>& qmatoms) {
            std::vector<std::string> UniqueElements= FindUniqueElements(qmatoms);
           
            BasisSet ecp;
            ecp.LoadPseudopotentialSet(_ecp_name);

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Loaded Pseudopotentials " << _ecp_name << flush;

            for (const std::string& element_name:UniqueElements) {
                       try{    
                        const Element& element = ecp.getElement(element_name);
                       }catch(std::runtime_error& error){
                         CTP_LOG(ctp::logDEBUG, *_pLog) << "No pseudopotential for " << element_name<<" available" << flush;
                         continue;
                       }
                       const Element& element = ecp.getElement(element_name);
                        // element name, [possibly indeces of centers], zero to indicate the end
                        com_file << element_name << " 0\n"
                                << _ecp_name << " "
                                << element.getLmax() << " " << element.getNcore() << endl;

                       for (const Shell& shell:element) {
                            // shell type, number primitives, scale factor
                            com_file << shell.getType() << endl;
                            com_file << shell.getSize() << endl;

                           for (const GaussianPrimitive& gaussian:shell) {
                                com_file << gaussian._power << " " << FortranFormat(gaussian._decay) << " " << FortranFormat(gaussian._contraction[0]) << endl;
                            }
                        }
                    }     
            com_file << endl;
            return;
        }

        /* For QM/MM the molecules in the MM environment are represented by
         * their atomic partial charge distributions. Triggered by the option
         * keyword "charge" Gaussian expects them in x,y,z,q format in the
         * input file. In g03 AFTER basis sets and ECPs, in g09 BEFORE.
         */
     
        void Gaussian::WriteBackgroundCharges(std::ofstream& com_file) {
            
            boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
            for (std::shared_ptr<ctp::PolarSeg> seg:_PolarSegments) {
                for (ctp::APolarSite* site:*seg) {
                    
                    string sitestring=boost::str(fmt % ((site->getPos().getX())*votca::tools::conv::nm2ang) 
                            % (site->getPos().getY()*votca::tools::conv::nm2ang) 
                            % (site->getPos().getZ()*votca::tools::conv::nm2ang) 
                            % site->getQ00());
                    if (site->getQ00() != 0.0) com_file << sitestring << endl;

                    if (site->getRank() > 0 || _with_polarization ) {

                        std::vector< std::vector<double> > _split_multipoles = SplitMultipoles(site);
                        for (const auto& mpoles:_split_multipoles){
                           string multipole=boost::str( fmt % mpoles[0] % mpoles[1] % mpoles[2] % mpoles[3]);
                            com_file << multipole << endl;

                        }
                    }
                }
            }
            com_file << endl;
            return;
        }

        /* An initial guess for the electron density can be provided by
         * a set of molecular orbital coefficients in the input file,
         * triggered by the 'guess=cards' keyword. This MUST be done in
         * Fortran fixed format 5D15.8. The information about the guess
         * itself is taken from a prepared orbitals object.
         */
        void Gaussian::WriteGuess(Orbitals& orbitals_guess, std::ofstream& com_file) {

            std::vector<int> sort_index = orbitals_guess.SortEnergies();
            ReorderMOsBack(orbitals_guess);

            com_file << "(5D15.8)" << endl;

            int level = 1;
            int ncolumns = 5;

            for (std::vector< int > ::iterator soi = sort_index.begin(); soi != sort_index.end(); ++soi) {


                com_file << setw(5) << level << endl;

                Eigen::VectorXd mr = orbitals_guess.MOCoefficients().col(*soi);

                int column = 1;
                for (unsigned j = 0; j < mr.size(); ++j) {
                    com_file << FortranFormat(mr[j]);
                    if (column == ncolumns) {
                        com_file << std::endl;
                        column = 0;
                    }
                    column++;
                }

                level++;
                if (column != 1) com_file << endl;
            }
            com_file << 0 << endl;

            return;
        }

        /* For output of the AO matrix of Vxc using the patched g03 version,
         * g03 has to be called a second time after completing the single-point
         * SCF calculation. A second input file is generated based on the
         * originally specified options by forcing to read the converged
         * electron density from the checkpoint file, setting run to serial.
         */
        void Gaussian::WriteVXCRunInputFile() {
            std::ofstream com_file2;

            std::string com_file_name_full2 = _run_dir + "/" + _input_vxc_file_name;

            com_file2.open(com_file_name_full2.c_str());
            // header
            if (_chk_file_name.size()) com_file2 << "%chk=" << _chk_file_name << endl;
            if (_memory.size()) com_file2 << "%mem=" << _memory << endl;
            com_file2 << "%nprocshared=1" << endl;

            // adjusting the options line to Vxc output only
            std::string options_vxc = _options;
            boost::algorithm::replace_all(options_vxc, "pseudo=read", "Geom=AllCheck");
            boost::algorithm::replace_all(options_vxc, "/gen", " chkbasis");
            boost::algorithm::replace_all(options_vxc, "punch=mo", "guess=read");
            boost::algorithm::replace_all(options_vxc, "guess=tcheck", "");
            boost::algorithm::replace_all(options_vxc, "guess=huckel", "");
            boost::algorithm::replace_all(options_vxc, "charge", "charge=check");
            if (options_vxc.size()) com_file2 << options_vxc << endl;

            com_file2 << endl;
            com_file2 << "VXC output run \n";
            com_file2 << endl;
            com_file2.close();
            return;
        }


        /* Coordinates are written in standard Element,x,y,z format to the
         * input file.
         */
        void Gaussian::WriteCoordinates(std::ofstream& com_file, std::vector<QMAtom*>& qmatoms) {
            std::vector< QMAtom* >::iterator it;
            for (it = qmatoms.begin(); it < qmatoms.end(); it++) {
              tools::vec pos=(*it)->getPos()*tools::conv::bohr2ang;
                    com_file << setw(3) << (*it)->getType().c_str()
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX()
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY()
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ()
                            << endl;
               
            }
            com_file << endl;
            return;
        }

        /* Standard Gaussian Header is written to the input file, with checkpoint,
         * memory, shared processor request, option string containing all
         * relevant keywords, charge, and spin information.
         */
        void Gaussian::WriteHeader(std::ofstream& com_file) {
            if (_chk_file_name.size()) com_file << "%chk=" << _chk_file_name << endl;
            if (_memory.size()) com_file << "%mem=" << _memory << endl;
            if (_threads > 0) com_file << "%nprocshared=" << _threads << endl;
            if (_options.size()) com_file << _options << endl;

            com_file << endl;
            com_file << "TITLE ";

            com_file << endl << endl;
            com_file << setw(2) << _charge << setw(2) << _spin << endl;
            return;
        }

        /**
         * Prepares the com file from a vector of segments
         * Appends a guess constructed from monomer orbitals if supplied
         */
        bool Gaussian::WriteInputFile(Orbitals& orbitals) {

            std::string temp_suffix = "/id";
            std::string scratch_dir_backup = _scratch_dir;

            std::ofstream com_file;
            std::string _com_file_name_full = _run_dir + "/" + _input_file_name;
            com_file.open(_com_file_name_full.c_str());

            // header
            WriteHeader(com_file);

            std::vector< QMAtom* > qmatoms = orbitals.QMAtoms();

            WriteCoordinates(com_file, qmatoms);

            /* The order in which the following information has to appear
             * in the Gaussian input file is different from version g03 to
             * version g09. The newest version (g2016) is not supported.
             */
            if (_executable == "g03") {

                // if we need to write basis sets, do it now
                if (_write_basis_set) WriteBasisset(com_file, qmatoms);

                // write ECPs
                if (_write_pseudopotentials) WriteECP(com_file, qmatoms);

                // write the background charges
                //if (_write_charges) WriteBackgroundCharges(_com_file, qmatoms);
                if (_write_charges) WriteBackgroundCharges(com_file);

                // write inital guess
                if (_write_guess){
                    WriteGuess(orbitals, com_file);
                }

            } else if (_executable == "g09") {

                // write the background charges
                //if (_write_charges) WriteBackgroundCharges(_com_file, qmatoms);
                if (_write_charges) WriteBackgroundCharges(com_file);
                // if we need to write basis sets, do it now
                if (_write_basis_set) WriteBasisset(com_file, qmatoms);

                // write ECPs
                if (_write_pseudopotentials) WriteECP(com_file, qmatoms);

                // write inital guess
                if (_write_guess){
                    WriteGuess(orbitals, com_file);
                }

            } else {
                throw std::runtime_error("Gaussian executable unknown. Must be either g03 or g09.");
            }

            // for Vxc AO matrix output only with pre-compiled G03
            if (_output_Vxc) WriteVXCRunInputFile();


            com_file << _vdWfooter << endl;

            com_file << endl;
            com_file.close();
            // and now generate a shell script to run both jobs
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

            _scratch_dir = scratch_dir_backup + temp_suffix;
            WriteShellScript();
            _scratch_dir = scratch_dir_backup;

            return true;
        }

        /* Gaussian will be executed within a shell in order to set some
         * environment variables for the local SCRATCH directory and
         * (legacy mode) running a second instance for AO matrix of Vxc
         * using patched g03. This function writes the shell script.
         */
        bool Gaussian::WriteShellScript() {
            std::ofstream shell_file;

            std::string shell_file_name_full = _run_dir + "/" + _shell_file_name;

            shell_file.open(shell_file_name_full.c_str());

            shell_file << "#!/bin/tcsh" << endl;
            shell_file << "mkdir -p " << _scratch_dir << endl;
            shell_file << "setenv GAUSS_SCRDIR " << _scratch_dir << endl;
            shell_file << _executable << " " << _input_file_name << endl;
            if (_output_Vxc) {
                shell_file << "rm fort.22" << endl;
                shell_file << "setenv DoPrtXC YES" << endl;
                shell_file << _executable << " " << _input_vxc_file_name << " >& /dev/null " << endl;
                shell_file << "setenv DoPrtXC NO" << endl;
                shell_file << "rm $GAUSS_SCRDIR/*" << endl;
            }
            shell_file.close();

            return true;
        }

        /**
         * Runs the Gaussian job.
         */
        bool Gaussian::Run( Orbitals& orbitals ) {
            CTP_LOG(ctp::logDEBUG, *_pLog) << "GAUSSIAN: running [" << _executable << " " << _input_file_name << "]" << flush;

            if (std::system(NULL)) {
                // if scratch is provided, run the shell script;
                // otherwise run gaussian directly and rely on global variables
                std::string command;
                if (_scratch_dir.size() != 0 || _output_Vxc) {
                    command = "cd " + _run_dir + "; tcsh " + _shell_file_name;
                    //            _command  = "cd " + _run_dir + "; mkdir -p " + _scratch_dir +"; " + _executable + " " + _input_file_name;
                } else {
                    command = "cd " + _run_dir + "; mkdir -p $GAUSS_SCRDIR; " + _executable + " " + _input_file_name;
                }
                int check = std::system(command.c_str());
                if (check == -1) {
                    CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                    return false;
                }
                if (CheckLogFile()) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "GAUSSIAN: finished job" << flush;
                    return true;
                } else {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "GAUSSIAN: job failed" << flush;
                }
            } else {
                CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                return false;
            }
            return true;
        }

        /**
         * Cleans up after the Gaussian job
         */
        void Gaussian::CleanUp() {

            // cleaning up the generated files
            if (_cleanup.size() != 0) {

                CTP_LOG(ctp::logDEBUG, *_pLog) << "Removing " << _cleanup << " files" << flush;
                tools::Tokenizer tok_cleanup(_cleanup, ", ");
                std::vector <std::string> cleanup_info;
                tok_cleanup.ToVector(cleanup_info);
                for (const std::string& substring:cleanup_info) {

                    if (substring == "com") {
                        std::string file_name = _run_dir + "/" + _input_file_name;
                        remove(file_name.c_str());
                        if (_output_Vxc) {
                            std::string file_name = _run_dir + "/" + _input_vxc_file_name;
                            remove(file_name.c_str());
                        }
                    }

                    if (substring == "sh") {
                        std::string file_name = _run_dir + "/" + _shell_file_name;
                        remove(file_name.c_str());
                    }

                    if (substring == "log") {
                        std::string file_name = _run_dir + "/" + _log_file_name;
                        remove(file_name.c_str());
                        if (_output_Vxc) {
                            size_t lastdot = _log_file_name.find_last_of(".");
                            if (lastdot == std::string::npos) {
                                cerr << endl;
                                cerr << "Could not remove Vxc log file" << flush;
                            }
                            std::string file_name2 = file_name.substr(0, lastdot) + "-2.log";
                            remove(file_name2.c_str());
                        }
                    }

                    if (substring == "chk") {
                        std::string file_name = _run_dir + "/" + _chk_file_name;
                        remove(file_name.c_str());
                    }

                    if (substring == "fort.7") {
                        std::string file_name = _run_dir + "/" + substring;
                        remove(file_name.c_str());
                        if (_output_Vxc) {
                            std::string file_name = _run_dir + "/" + "fort.24";
                            remove(file_name.c_str());
                        }
                    }

                    if (substring == "gbs" && _write_basis_set) {
                        std::vector<std::string> fileswithfileending;
                        boost::filesystem::recursive_directory_iterator fit(_run_dir);
                        boost::filesystem::recursive_directory_iterator endit;
                        while (fit != endit) {
                            if (boost::filesystem::is_regular_file(* fit) && fit->path().extension() == substring) fileswithfileending.push_back(fit->path().filename().string());
                            ++fit;
                        }
                        for (const auto filename : fileswithfileending) {
                            std::string file_name = _run_dir + "/" + filename;
                            remove(file_name.c_str());
                        }
                    }

                }
            }
            return;
        }

        /**
         * Reads in the MO coefficients from a GAUSSIAN fort.7 file
         */
        bool Gaussian::ParseOrbitalsFile(Orbitals & orbitals) {
            std::map <int, std::vector<double> > coefficients;
            std::map <int, double> energies;

            std::string line;
            unsigned levels = 0;
            unsigned level = 0;
            unsigned basis_size = 0;

            std::string orb_file_name_full = _orb_file_name;
            if (_run_dir != "") orb_file_name_full = _run_dir + "/" + _orb_file_name;
            std::ifstream input_file(orb_file_name_full.c_str());

            if (input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "File " << _orb_file_name << " with molecular orbitals is not found " << flush;
                return false;
            } else {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Reading MOs from " << _orb_file_name << flush;
            }

            // number of coefficients per line is  in the first line of the file (5D15.8)
            getline(input_file, line);
            std::vector<std::string> strs;
            boost::algorithm::split(strs, line, boost::is_any_of("(D)"));
            std::string format = strs.at(2);


            while (input_file) {

                getline(input_file, line);
                // if a line has an equality sign, must be energy
                std::string::size_type energy_pos = line.find("=");

                if (energy_pos != std::string::npos) {

                    std::vector<std::string> results;
                    boost::trim(line);

                    boost::algorithm::split(results, line, boost::is_any_of("\t ="),
                            boost::algorithm::token_compress_on);
                    //cout << results[1] << ":" << results[2] << ":" << results[3] << ":" << results[4] << endl;

                    level = boost::lexical_cast<int>(results.front());
                    boost::replace_first(results.back(), "D", "e");
                    energies[ level ] = boost::lexical_cast<double>(results.back());
                    levels++;

                } else {

                    while (line.size() > 1) {
                        std::string coefficient;
                        coefficient.assign(line, 0, 15);
                        boost::trim(coefficient);
                        boost::replace_first(coefficient, "D", "e");
                        double coef = boost::lexical_cast<double>(coefficient);
                        coefficients[ level ].push_back(coef);
                        line.erase(0, 15);
                    }
                }
            }

            // some sanity checks
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Energy levels: " << levels << flush;

            std::map< int, std::vector<double> >::iterator iter = coefficients.begin();
            basis_size = iter->second.size();

            for (iter = coefficients.begin()++; iter != coefficients.end(); iter++) {
                if (iter->second.size() != basis_size) {
                    CTP_LOG(ctp::logERROR, *_pLog) << "Error reading " << _orb_file_name << ". Basis set size change from level to level." << flush;
                    return false;
                }
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Basis set size: " << basis_size << flush;

            // copying information to the orbitals object
            orbitals.setBasisSetSize(basis_size); // = _basis_size;

            // copying energies to the orbitals object
           Eigen::VectorXd &mo_energies = orbitals.MOEnergies();
            mo_energies.resize(levels);
            for (int i = 0; i < mo_energies.size(); i++) mo_energies[i] = energies[ i + 1 ];

            // copying mo coefficients to the orbitals object
            Eigen::MatrixXd &mo_coefficients = orbitals.MOCoefficients();
            mo_coefficients.resize(levels, basis_size);
            for (int i = 0; i < mo_coefficients.rows(); i++){
                for (int j = 0; j < mo_coefficients.cols(); j++){
                    mo_coefficients(j, i) = coefficients[i + 1][j];
                }
            }
            
            ReorderOutput(orbitals);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "GAUSSIAN: done reading MOs" << flush;

            return true;
        }

        bool Gaussian::CheckLogFile() {

            // check if the log file exists
            boost::filesystem::path arg_path;
            char ch;

            std::string full_name = (arg_path / _run_dir / _log_file_name).c_str();
            ifstream input_file(full_name.c_str());

            if (input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "GAUSSIAN: " << full_name << " is not found" << flush;
                return false;
            };

            input_file.seekg(0, ios_base::end); // go to the EOF

            // get empty lines and end of lines out of the way
            do {
                input_file.seekg(-2, ios_base::cur);
                input_file.get(ch);
            } while (ch == '\n' || ch == ' ' || ch == '\t' || (int) input_file.tellg() == -1);

            // get the beginning of the line or the file
            do {
                input_file.seekg(-2, ios_base::cur);
                input_file.get(ch);
            } while (ch != '\n' && (int) input_file.tellg() != -1);

            std::string line;
            getline(input_file, line); 
            input_file.close();

            std::string::size_type self_energy_pos = line.find("Normal termination of Gaussian");
            if (self_energy_pos == std::string::npos) {
                CTP_LOG(ctp::logERROR, *_pLog) << "GAUSSIAN: " << full_name << " is incomplete" << flush;
                return false;
            } else {
                return true;
            }
        }

        bool Gaussian::ReadESPCharges(Orbitals& orbitals, std::string& line, ifstream& input_file){
          std::string::size_type charge_pos = line.find("Charges from ESP fit, RMS");
          bool has_charges=false;
          if (charge_pos != std::string::npos && _get_charges) {
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting charges" << flush;
            has_charges = true;
            getline(input_file, line);
            getline(input_file, line);
            
            bool _has_atoms = orbitals.hasQMAtoms();
            
            std::vector<std::string> row=GetLineAndSplit(input_file, "\t ");
            int nfields = row.size();
            
            while (nfields == 3) {
              int atom_id = boost::lexical_cast< int >(row.at(0));
              std::string atom_type = row.at(1);
              double atom_charge = boost::lexical_cast< double >(row.at(2));
              row=GetLineAndSplit(input_file, "\t ");
              nfields = row.size();
              QMAtom* pAtom;
              if (_has_atoms == false) {
                pAtom =orbitals.AddAtom(atom_id - 1,atom_type, 0, 0, 0);
              } else {
                pAtom = orbitals.QMAtoms().at(atom_id - 1);
              }
              pAtom->setPartialcharge(atom_charge);
            }
          }
          return has_charges;
        }

        /**
         * Parses the Gaussian Log file and stores data in the Orbitals object
         */
        bool Gaussian::ParseLogFile(Orbitals & orbitals) {
            std::string line;
            std::vector<std::string> results;
            bool has_occupied_levels = false;
            bool has_unoccupied_levels = false;
            bool has_number_of_electrons = false;
            bool has_basis_set_size = false;
            bool has_overlap_matrix = false;
            bool has_charges = false;
            bool has_self_energy = false;

            bool read_vxc = false;

            int occupied_levels = 0;
            int unoccupied_levels = 0;
            int number_of_electrons = 0;
            int basis_set_size = 0;
            int cart_basis_set_size = 0;

            CTP_LOG(ctp::logDEBUG, *_pLog) << "GAUSSIAN: parsing " << _log_file_name << flush;

            std::string log_file_name_full = _log_file_name;
            if (_run_dir != "") log_file_name_full = _run_dir + "/" + _log_file_name;

            // check if LOG file is complete
            if (!CheckLogFile()) return false;

            // save qmpackage name
            orbitals.setQMpackage("gaussian");
            orbitals.setDFTbasis(_basisset_name);


            if (_write_pseudopotentials) {
                orbitals.setECP(_ecp_name);
            }

            read_vxc = _output_Vxc;
            bool vxc_found = false;
            // Start parsing the file line by line
            ifstream input_file(log_file_name_full.c_str());
            while (input_file) {

                getline(input_file, line);
                boost::trim(line);

                /* Check for ScaHFX = factor of HF exchange included in functional */
                std::string::size_type HFX_pos = line.find("ScaHFX=");
                if (HFX_pos != std::string::npos) {
                    boost::algorithm::split(results, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    double ScaHFX = boost::lexical_cast<double>(results.back());
                    orbitals.setScaHFX(ScaHFX);
                    vxc_found = true;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "DFT with " << ScaHFX << " of HF exchange!" << flush;
                }

                /*
                 * number of occupied and virtual orbitals
                 * N alpha electrons      M beta electrons
                 */
                std::string::size_type electrons_pos = line.find("alpha electrons");
                if (electrons_pos != std::string::npos) {
                    boost::algorithm::split(results, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    has_number_of_electrons = true;
                    number_of_electrons = boost::lexical_cast<int>(results.front());
                    orbitals.setNumberOfElectrons(number_of_electrons);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Alpha electrons: " << number_of_electrons << flush;
                }

                /*
                 * basis set size
                 * N basis functions,  M primitive gaussians,   K cartesian basis functions
                 */
                std::string::size_type basis_pos = line.find("basis functions,");
                if (basis_pos != std::string::npos) {
                    boost::algorithm::split(results, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    has_basis_set_size = true;
                    basis_set_size = boost::lexical_cast<int>(results.front());
                    orbitals.setBasisSetSize(basis_set_size);
                    cart_basis_set_size = boost::lexical_cast<int>(results[6]);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Basis functions: " << basis_set_size << flush;
                    if (read_vxc) {
                        CTP_LOG(ctp::logDEBUG, *_pLog) << "Cartesian functions: " << cart_basis_set_size << flush;
                    }
                }

                /*
                 * energies of occupied/unoccupied levels
                 * Alpha  occ.(virt.) eigenvalues -- e1 e2 e3 e4 e5
                 */
                std::string::size_type eigenvalues_pos = line.find("Alpha");
                if (eigenvalues_pos != std::string::npos) {

                    std::list<std::string> stringList;

                    while (eigenvalues_pos != std::string::npos && !has_occupied_levels && !has_unoccupied_levels) {

                        boost::iter_split(stringList, line, boost::first_finder("--"));

                        std::vector<std::string> energies;
                        boost::trim(stringList.back());

                        boost::algorithm::split(energies, stringList.back(), boost::is_any_of("\t "), boost::algorithm::token_compress_on);

                        if (stringList.front().find("virt.") != std::string::npos) {
                            unoccupied_levels += energies.size();
                            energies.clear();
                        }

                        if (stringList.front().find("occ.") != std::string::npos) {
                            occupied_levels += energies.size();
                            energies.clear();
                        }

                        getline(input_file, line);
                        eigenvalues_pos = line.find("Alpha");
                        boost::trim(line);

                        if (eigenvalues_pos == std::string::npos) {
                            has_occupied_levels = true;
                            has_unoccupied_levels = true;
                            orbitals.setNumberOfLevels(occupied_levels, unoccupied_levels);
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "Occupied levels: " << occupied_levels << flush;
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "Unoccupied levels: " << unoccupied_levels << flush;
                        }
                    } // end of the while loop
                } // end of the eigenvalue parsing

                /*
                 *  Partial charges from the input file
                 */
                has_charges=ReadESPCharges(orbitals, line,input_file);

                /*
                 * Coordinates of the final configuration
                 * stored in the archive at the end of the file
                 */
                int cpn = 0; // marker appearence marker
                std::string::size_type coordinates_pos = line.find("\\");

                if (coordinates_pos != std::string::npos && cpn == 0) {
                    ++cpn; // updates but ignores
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting the coordinates" << flush;
                    boost::trim(line);
                    std::string archive = line;
                    while (line.size() != 0) {
                        getline(input_file, line);
                        boost::trim(line);
                        archive += line;
                    }

                    bool has_atoms = orbitals.hasQMAtoms();
                    std::list<std::string> stringList;
                    std::vector<std::string> results;
                    boost::iter_split(stringList, archive, boost::first_finder("\\\\"));

                    std::list<std::string>::iterator coord_block = stringList.begin();
                    std::advance(coord_block, 3);

                    std::vector<std::string> atom_block;
                    boost::algorithm::split(atom_block, *coord_block, boost::is_any_of("\\"), boost::algorithm::token_compress_on);

                    std::vector<std::string>::iterator atom_block_it;
                    int aindex = 0;

                    for (atom_block_it = ++atom_block.begin(); atom_block_it != atom_block.end(); ++atom_block_it) {
                        std::vector<std::string> atom;
                        boost::algorithm::split(atom, *atom_block_it, boost::is_any_of(","), boost::algorithm::token_compress_on);
                        std::string atom_type = atom.front();
                        std::vector<std::string>::iterator it_atom;
                        it_atom = atom.end();
                        double z = boost::lexical_cast<double>(*(--it_atom));
                        double y = boost::lexical_cast<double>(*(--it_atom));
                        double x = boost::lexical_cast<double>(*(--it_atom));
                        tools::vec pos=tools::vec(x,y,z);
                        pos*=tools::conv::ang2bohr;

                        if (has_atoms == false) {
                            orbitals.AddAtom(aindex,atom_type, pos);
                        } else {
                            QMAtom* pAtom = orbitals.QMAtoms().at(aindex);
                            pAtom->setPos(pos);
                            
                        }
                        aindex++;
                    }
                    // get the QM energy out
                    std::advance(coord_block, 1);
                    std::vector<std::string> block;
                    std::vector<std::string> energy;
                    boost::algorithm::split(block, *coord_block, boost::is_any_of("\\"), boost::algorithm::token_compress_on);
                    map<std::string, std::string> properties;
                    std::vector<std::string>::iterator block_it;
                    for (block_it = block.begin(); block_it != block.end(); ++block_it) {
                        std::vector<std::string> property;
                        boost::algorithm::split(property, *block_it, boost::is_any_of("="), boost::algorithm::token_compress_on);
                        properties[property[0]] = property[1];
                    }
                    if (properties.count("HF") > 0) {
                        double energy_hartree = boost::lexical_cast<double>(properties["HF"]);
                        orbitals. setQMEnergy(tools::conv::hrt2ev * energy_hartree);
                        CTP_LOG(ctp::logDEBUG, *_pLog) << (boost::format("QM energy[eV]: %4.6f ") % orbitals.getQMEnergy()).str() << flush;
                    } else {
                        cout << endl;
                        throw std::runtime_error("ERROR No energy in archive");
                    }

                }

                std::string::size_type self_energy_pos = line.find("Self energy of the charges");

                if (self_energy_pos != std::string::npos) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting the self energy\n";
                    std::vector<std::string> block;
                    std::vector<std::string> energy;
                    boost::algorithm::split(block, line, boost::is_any_of("="), boost::algorithm::token_compress_on);
                    boost::algorithm::split(energy, block[1], boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    orbitals.setSelfEnergy(tools::conv::hrt2ev * boost::lexical_cast<double> (energy[1]));
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Self energy " << orbitals.getSelfEnergy() << flush;

                }
 
                std::string::size_type overlap_pos = line.find("*** Overlap ***");
                if (overlap_pos != std::string::npos) {

                    // prepare the container
                    Eigen::MatrixXd& overlap=orbitals.AOOverlap();
                    overlap.resize(basis_set_size,basis_set_size);
                    has_overlap_matrix = true;
                    std::vector<int> j_indeces;
                    int n_blocks = 1 + ((basis_set_size - 1) / 5);
                    getline(input_file, line);
                    boost::trim(line);

                    for (int _block = 0; _block < n_blocks; _block++) {
                        // first line gives the j index in the matrix
                        boost::tokenizer<> tok(line);
                        std::transform(tok.begin(), tok.end(), std::back_inserter(j_indeces), &boost::lexical_cast<int, std::string>);

                        // read the block of max _basis_size lines + the following header
                        for (int i = 0; i <= basis_set_size; i++) {
                            getline(input_file, line);
                            if (std::string::npos == line.find("D")) break;
                            // split the line on the i index and the rest
                            std::vector<std::string> row=GetLineAndSplit(input_file, "\t ");
                            int i_index = boost::lexical_cast<int>(row.front());
                            row.erase(row.begin());
                            std::vector<int>::iterator j_iter = j_indeces.begin();

                            for (std::string& coefficient: row) {
                                boost::replace_first(coefficient, "D", "e");
                                int j_index = *j_iter;
                                overlap(i_index - 1, j_index - 1) = boost::lexical_cast<double>(coefficient);
                                overlap(j_index - 1, i_index - 1) = boost::lexical_cast<double>(coefficient);
                                j_iter++;
                            }
                        }
                        // clear the index for the next block
                        j_indeces.clear();
                    } // end of the blocks
 
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Read the overlap matrix" << flush;
                } // end of the if "Overlap" found
                // check if all information has been accumulated and quit
                if (has_number_of_electrons &&
                        has_basis_set_size &&
                        has_occupied_levels &&
                        has_unoccupied_levels &&
                        has_overlap_matrix &&
                        has_charges &&
                        has_self_energy
                        ) break;

            } // end of reading the file line-by-line

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
            input_file.close();

            if (!vxc_found) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "WARNING === WARNING \n, could not find ScaHFX= entry in log."
                        "\n probably you forgt #P in the beginning of the input file.\n"
                        " If you are running a hybrid functional calculation redo it! Now! Please!\n ===WARNING=== \n"
                        << flush;
                orbitals.setScaHFX(0.0);
            }
            // - parse atomic orbitals Vxc matrix
            if (read_vxc) {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Parsing fort.24 for Vxc" << flush;
                std::string log_file_name_full;
                if (_run_dir == "") {
                    log_file_name_full = "fort.24";
                } else {
                    log_file_name_full = _run_dir + "/fort.24";
                }

                ifstream input_file(log_file_name_full.c_str());
                if (input_file.good()) {
                    // prepare the container
                    Eigen::MatrixXd vxc=Eigen::MatrixXd::Zero(cart_basis_set_size,cart_basis_set_size);
                    std::vector<int> j_indeces;
                    // Start parsing the file line by line

                    while (input_file) {
                        getline(input_file, line);
                        if (input_file.eof()) break;

                        std::vector<std::string> row;
                        boost::trim(line);
                        boost::algorithm::split(row, line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);

                        int i_index = boost::lexical_cast<int>(row[0]);
                        int j_index = boost::lexical_cast<int>(row[1]);
                        vxc(i_index - 1, j_index - 1) = boost::lexical_cast<double>(row[2]);
                        vxc(j_index - 1, i_index - 1) = boost::lexical_cast<double>(row[2]);
                    }

                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
                    input_file.close();
                BasisSet dftbasisset;
                dftbasisset.LoadBasisSet(_basisset_name);
                if(!orbitals.hasQMAtoms()){
                    throw runtime_error("Orbitals object has no QMAtoms");
                }
                AOBasis dftbasis;
                dftbasis.AOBasisFill(dftbasisset, orbitals.QMAtoms());
                Eigen::MatrixXd carttrafo=dftbasis.getTransformationCartToSpherical(getPackageName());
                orbitals.AOVxc()=carttrafo*vxc*carttrafo.transpose();
                } else {
                    throw std::runtime_error("Vxc file does not exist.");
                }
            }
            return true;
        }

        std::string Gaussian::FortranFormat(double number) {
            std::stringstream ssnumber;
            if (number >= 0) ssnumber << " ";
            ssnumber << setiosflags(ios::fixed) << setprecision(8) << std::scientific << number;
            std::string snumber = ssnumber.str();
            boost::replace_first(snumber, "e", "D");
            return snumber;
        }


    }
}
