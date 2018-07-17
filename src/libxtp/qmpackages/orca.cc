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

#include "orca.h"
#include <votca/xtp/qminterface.h>

#include <votca/ctp/segment.h>

#include <votca/tools/elements.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <stdio.h>
#include <iomanip>
#include <sys/stat.h>
#include <vector>



namespace votca {
    namespace xtp {
      using namespace std;

        void Orca::Initialize(tools::Property *options) {

            //good luck

            // Orca file names
            std::string fileName = "system";

            _xyz_file_name = fileName + ".xyz";
            _input_file_name = fileName + ".inp";
            _log_file_name = fileName + ".log";
            _shell_file_name = fileName + ".sh";
            _orb_file_name = fileName + ".gbw";

            std::string key = "package";
            std::string _name = options->get(key + ".name").as<std::string> ();

            if (_name != "orca") {
                cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error("Wrong options file");
            }

            _executable = options->get(key + ".executable").as<std::string> ();
            _charge = options->get(key + ".charge").as<int> ();
            _spin = options->get(key + ".spin").as<int> ();
            _options = options->get(key + ".options").as<std::string> ();
            _memory = options->get(key + ".memory").as<std::string> ();
            _threads = options->get(key + ".threads").as<int> ();
            _scratch_dir = options->get(key + ".scratch").as<std::string> ();
            _cleanup = options->get(key + ".cleanup").as<std::string> ();
            _auxbasisset_name = options->get(key + ".auxbasisset").as<std::string> ();



            if (options->exists(key + ".outputVxc")) {
                _output_Vxc = options->get(key + ".outputVxc").as<bool> ();
            } else _output_Vxc = false;
            if (_output_Vxc) {
                throw std::runtime_error("Sorry " + _name + " does not support Vxc output");
            }

            _basisset_name = options->get(key + ".basisset").as<std::string> ();
            _write_basis_set = options->get(key + ".writebasisset").as<bool> ();
            _write_pseudopotentials = options->get(key + ".writepseudopotentials").as<bool> ();
            if ( _write_pseudopotentials )  _ecp_name = options->get(key + ".ecp").as<std::string> ();


            // check if the optimize keyword is present, if yes, read updated coords
            std::string::size_type iop_pos = _options.find(" Opt"); /*optimization word in orca*/
            if (iop_pos != std::string::npos) {
                _is_optimization = true;
            } else {
                _is_optimization = false;
            }

            // check if the esp keyword is present, if yes, get the charges and save them
            iop_pos = _options.find(" chelpg"); /*electrostatic potential related to partial atomic charges I guess is chelpg in orca but check */
            if (iop_pos != std::string::npos ||  _options.find(" CHELPG")!= std::string::npos) {
                _get_charges = true;
            } else {
                _get_charges = false;
            }


            // check if the guess should be prepared, if yes, append the guess later
            _write_guess = false;
            iop_pos = _options.find("Guess MORead");
            if (iop_pos != std::string::npos) _write_guess = true;
            iop_pos = _options.find("Guess MORead\n");
            if (iop_pos != std::string::npos) _write_guess = true;
        }


        /* Custom basis sets are written on a per-element basis to
         * the system.bas/aux file(s), which are then included in the
         * Orca input file using GTOName = "system.bas/aux"
         */
        void Orca::WriteBasisset(std::vector<QMAtom*>& qmatoms, std::string& _bs_name, std::string& _el_file_name) {

            tools::Elements _elements;
            list<std::string> elements;
            BasisSet bs;
            bs.LoadBasisSet(_bs_name);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Loaded Basis Set " << _bs_name << flush;
            ofstream _el_file;

            _el_file.open(_el_file_name.c_str());
            //_com_file << "@" << "system.bas" << endl;
            _el_file << "$DATA" << endl;
            std::vector< QMAtom* >::iterator it;

            for (it = qmatoms.begin(); it < qmatoms.end(); it++) {
                std::string element_name = (*it)->getType();
                list<std::string>::iterator ite;
                ite = find(elements.begin(), elements.end(), element_name);
                if (ite == elements.end()) {
                    elements.push_back(element_name);
                    Element* element = bs.getElement(element_name);
                    _el_file << _elements.getEleFull(element_name) << endl;
                    for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                        Shell* shell = (*its);

                        string type = shell->getType();
                        // check combined shells

                        for (unsigned i = 0; i < type.size(); ++i) {
                            string subtype = string(type, i, 1);
                            _el_file << subtype << " " << shell->getSize() << endl;
                            int _sh_idx = 0;
                            for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                GaussianPrimitive* gaussian = *itg;
                                _sh_idx++;
                                _el_file << " " << _sh_idx << " " << indent(gaussian->decay);
                                _el_file << " " << indent(gaussian->contraction[FindLmax(subtype)]);

                                _el_file << endl;
                            }

                        }
                    }

                }
            }

            _el_file << "STOP\n";
            _el_file.close();

            return;
        }

        /* Coordinates are written in standard Element,x,y,z format to the
         * input file.
         */
        void Orca::WriteCoordinates(std::ofstream& _com_file, std::vector<QMAtom*>& qmatoms) {

            std::vector< QMAtom* >::iterator it;
            for (it = qmatoms.begin(); it < qmatoms.end(); it++) {
                tools::vec pos=(*it)->getPos()*tools::conv::bohr2ang;
                    _com_file << setw(3) << (*it)->getType().c_str()
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX()
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY()
                            << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ()
                            << endl;
            }
            _com_file << "* \n" << endl;
            return;
        }


        /* If custom ECPs are used, they need to be specified in the input file
         * in a section following the basis set includes.
         */
        void Orca::WriteECP(std::ofstream& _com_file, std::vector<QMAtom*>& qmatoms){

            std::vector< QMAtom* >::iterator it;


            _com_file << endl;

            list<std::string> elements;
            elements.push_back("H");
            elements.push_back("He");
            BasisSet ecp;
            ecp.LoadPseudopotentialSet(_ecp_name);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Loaded Pseudopotentials " <<_ecp_name << flush;



            for (it = qmatoms.begin(); it < qmatoms.end(); it++) {
                
                    std::string element_name = (*it)->getType();


                    list<std::string>::iterator ite;
                    ite = find(elements.begin(), elements.end(), element_name);
                    if (ite == elements.end()) {
                        elements.push_back(element_name);
                        Element* element = ecp.getElement(element_name);
                        _com_file << "\n" << "NewECP" << " " << element_name << endl;
                        _com_file << "N_core" << " " << element->getNcore() << endl;
                        //lmaxnum2lmaxname
                        _com_file << "lmax" << " " << getLName(element->getLmax()) << endl;

                        //For Orca the order doesn't matter but let's write it in ascending order
                        // write remaining shells in ascending order s,p,d...
                        for (int i = 0; i <= element->getLmax(); i++) {
                            for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                Shell* shell = (*its);
                                if (shell->getLmax() == i) {
                                    // shell type, number primitives, scale factor
                                    _com_file << shell->getType() << " " << shell->getSize() << endl;
                                    int _sh_idx = 0;
                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                        _sh_idx++;
                                        _com_file << _sh_idx << " " << gaussian->decay << " " << gaussian->contraction[0] << " " << gaussian->power << endl;
                                    }
                                }
                            }
                        }
                        _com_file << "end\n " << "\n" << endl;
                    }

                
            }
            return;
        }

        /* For QM/MM the molecules in the MM environment are represented by
         * their atomic partial charge distributions. ORCA expects them in
         * q,x,y,z format in a separate file "background.crg"
         */
        void Orca::WriteBackgroundCharges() {
            std::ofstream _crg_file;
            std::string _crg_file_name_full = _run_dir + "/background.crg";
            _crg_file.open(_crg_file_name_full.c_str());
            int _total_background = 0;

            std::vector< ctp::PolarSeg* >::iterator it;
            for (it = _PolarSegments.begin(); it < _PolarSegments.end(); it++) {
                vector<ctp::APolarSite*> ::iterator pit;
                for (pit = (*it)->begin(); pit < (*it)->end(); ++pit) {
                    if ((*pit)->getQ00() != 0.0) _total_background++;

                    if ((*pit)->getRank() > 0 || _with_polarization ) {

                        std::vector<std::vector<double>> _split_multipoles = SplitMultipoles(*pit);
                        _total_background+= _split_multipoles.size();
                    }
                }
            } //counting only

            
            _crg_file << _total_background << endl;
            boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
            //now write
            for (it = _PolarSegments.begin(); it < _PolarSegments.end(); it++) {
                vector<ctp::APolarSite*> ::iterator pit;
                for (pit = (*it)->begin(); pit < (*it)->end(); ++pit) {
                    
                    string site=boost::str(fmt % (*pit)->getQ00() % (((*pit)->getPos().getX())*votca::tools::conv::nm2ang) 
                            % ((*pit)->getPos().getY()*votca::tools::conv::nm2ang) 
                            % ((*pit)->getPos().getZ()*votca::tools::conv::nm2ang) 
                            );
                    if ((*pit)->getQ00() != 0.0) _crg_file << site << endl;

                    if ((*pit)->getRank() > 0 || _with_polarization ) {

                        std::vector< std::vector<double> > _split_multipoles = SplitMultipoles(*pit);
                        for (const auto& mpoles:_split_multipoles){
                           string multipole=boost::str( fmt % mpoles[3] % mpoles[0] % mpoles[1] % mpoles[2] );
                            _crg_file << multipole << endl;

                        }
                    }
                }
            }
            
            return;
        }

        /**
         * Prepares the *.inp file from a vector of segments
         * Appends a guess constructed from monomer orbitals if supplied, Not implemented yet
         */
        bool Orca::WriteInputFile(std::vector< ctp::Segment* > segments, Orbitals* orbitals_guess ) {

            std::vector<std::string> results;
            std::string temp_suffix = "/id";
            std::string scratch_dir_backup = _scratch_dir;
            std::ofstream _com_file;

            std::string _com_file_name_full = _run_dir + "/" + _input_file_name;

            _com_file.open(_com_file_name_full.c_str());
            // header
            _com_file << "* xyz  " << _charge << " " << _spin << endl;

            std::vector< QMAtom* > qmatoms;
            // This is needed for the QM/MM scheme, since only orbitals have
            // updated positions of the QM region, hence vector<Segments*> is
            // NULL in the QMMachine and the QM region is also printed here
            if (_write_charges) {
                qmatoms = orbitals_guess->QMAtoms();
            } else {
                QMInterface qmmface;
                qmatoms = qmmface.Convert(segments);
            }

            // put coordinates
            WriteCoordinates(_com_file, qmatoms);

            // add parallelization info
            _com_file << "%pal\n " << "nprocs " << _threads << "\nend" << "\n" << endl;

            // basis set info
            if (_write_basis_set) {

                std::string _el_file_name = _run_dir + "/" + "system.bas";
                WriteBasisset(qmatoms, _basisset_name, _el_file_name);
                _com_file << "%basis\n " << endl;
                _com_file << "GTOName" << " " << "=" << "\"system.bas\";" << endl;

                if (_auxbasisset_name != "") {
                    std::string _aux_file_name = _run_dir + "/" + "system.aux";
                    WriteBasisset(qmatoms, _auxbasisset_name, _aux_file_name);
                    _com_file << "GTOAuxName" << " " << "=" << "\"system.aux\";" << endl;
                }
            } // write_basis set

            // ECPs
            /* WRITING ECP INTO system.inp FILE for ORCA**/
            if (_write_pseudopotentials) {
                WriteECP(_com_file, qmatoms);
            } // write pseudopotentials


            /* END   OF WRITING BASISSET/ECP INTO system.inp FILE for ORCA*************/
            _com_file << "end\n " << "\n" << endl; //This end is for the basis set block


            if (_write_charges) {
                WriteBackgroundCharges();
            }

            _com_file << _options << "\n";

           

            _com_file << endl;
            _com_file.close();

            // and now generate a shell script to run both jobs, if neccessary
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

            _scratch_dir = scratch_dir_backup + temp_suffix;
            WriteShellScript();
            _scratch_dir = scratch_dir_backup;

            return true;
        }
        
        
   

        bool Orca::WriteShellScript() {
            ofstream _shell_file;

            std::string _shell_file_name_full = _run_dir + "/" + _shell_file_name;

            _shell_file.open(_shell_file_name_full.c_str());

            _shell_file << "#!/bin/bash" << endl;
            _shell_file << "mkdir -p " << _scratch_dir << endl;
            
            if(_write_guess){
              if(!(boost::filesystem::exists( _run_dir + "/molA.gbw" ) && boost::filesystem::exists( _run_dir + "/molB.gbw" )  )){
              throw runtime_error("Using guess relies on a molA.gbw and a molB.gbw file being in the directory.");
            }
              
              _shell_file<<_executable<<"_mergefrag molA.gbw molB.gbw dimer.gbw > merge.log"<<endl;
            }

            if (_threads == 1) {
                _shell_file << _executable << " " << _input_file_name << " > " << _log_file_name << endl; //" 2> run.error" << endl;
            } else {
                _shell_file << _executable << " " << _input_file_name << " > " << _log_file_name << endl; // " 2> run.error" << endl;
            }
            _shell_file.close();

            return true;
        }

        /**
         * Runs the Orca job.
         */
        bool Orca::Run( Orbitals* _orbitals ) {

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Running Orca job" << flush;

            if (std::system(NULL)) {

                // Orca overrides input information, if *.db and *.movecs files are present
                // better trash the old version
                std::string file_name = _run_dir + "/system.db";
                remove(file_name.c_str());
                file_name = _run_dir + "/" + _log_file_name;
                remove(file_name.c_str());
                file_name = _run_dir + "/" + _orb_file_name;
                //remove ( file_name.c_str() );

                std::string _command;
                if (_threads == 1) {
                    _command = "cd " + _run_dir + "; sh " + _shell_file_name;
                } else {
                    _command = "cd " + _run_dir + "; sh " + _shell_file_name;
                }
                //CTP_LOG(logDEBUG,*_pLog) << _command << flush;
                int check = std::system(_command.c_str());
                if (check == -1) {
                    CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                    return false;
                }

                if (CheckLogFile()) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Finished Orca job" << flush;
                    return true;
                } else {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Orca job failed" << flush;
                }
            } else {
                CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                return false;
            }

            return true;
        }

        /**
         * Cleans up after the Orca job
         */
        void Orca::CleanUp() {

          if(_write_guess){
            remove((_run_dir + "/" +"molA.gbw").c_str());
            remove((_run_dir + "/" +"molB.gbw").c_str());
            remove((_run_dir + "/" +"dimer.gbw").c_str());
            
          }
          
            // cleaning up the generated files
            if (_cleanup.size() != 0) {
                tools::Tokenizer tok_cleanup(_cleanup, ",");
                std::vector <std::string> _cleanup_info;
                tok_cleanup.ToVector(_cleanup_info);

                std::vector<std::string> ::iterator it;
                
                

                for (it = _cleanup_info.begin(); it != _cleanup_info.end(); ++it) {
                    if (*it == "inp") {
                        std::string file_name = _run_dir + "/" + _input_file_name;
                        remove(file_name.c_str());
                    }

                    if (*it == "bas") {
                        std::string file_name = _run_dir + "/system.bas";
                        remove(file_name.c_str());
                    }

                    if (*it == "log") {
                        std::string file_name = _run_dir + "/" + _log_file_name;
                        remove(file_name.c_str());
                    }

                    if (*it == "gbw") {
                        std::string file_name = _run_dir + "/" + _orb_file_name;
                        remove(file_name.c_str());
                    }

                    if (*it == "ges") {
                        std::string file_name = _run_dir + "/system.ges";
                        remove(file_name.c_str());
                    }
                    if (*it == "prop") {
                        std::string file_name = _run_dir + "/system.prop";
                        remove(file_name.c_str());
                    }
                }
            }
            return;
        }

        bool Orca::ParseLogFile(Orbitals* _orbitals) {
            const double _conv_Hrt_eV = tools::conv::hrt2ev;

            _orbitals->setQMpackage("orca");
            _orbitals->setDFTbasis(_basisset_name);

            if (_write_pseudopotentials) {
                _orbitals->setECP(_ecp_name);
            } 

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Parsing " << _log_file_name << flush;
            // return true;
            std::string _log_file_name_full = _run_dir + "/" + _log_file_name;
            // check if LOG file is complete
            if (!CheckLogFile()) return false;

            //std::map <int, std::std::vector<double> > _coefficients;
            std::map <int, double> _energies;
            std::map <int, double> _occ;

            std::string _line;
            unsigned _levels = 0;
            //unsigned _level;
            //unsigned _basis_size = 0;
            int _number_of_electrons = 0;
            //bool _has_basis_dim = false;

            std::vector<std::string> results;

            std::ifstream _input_file(_log_file_name_full.c_str());

            if (_input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "File " << _log_file_name_full << " not found " << flush;
                return false;
            } else {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Reading Coordinates and occupationnumbers and energies from " << _log_file_name_full << flush;
            }


            //Coordinates of the final configuration depending on whether it is an optimization or not



            while (_input_file) {
                getline(_input_file, _line);
                boost::trim(_line);



                if (_is_optimization) {
                    throw runtime_error("Not implemented yet!");
                }
                bool _found_optimization = true;

                std::string::size_type coordinates_pos = _line.find("CARTESIAN COORDINATES (ANGSTROEM)");

                if (_found_optimization && coordinates_pos != std::string::npos) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting the coordinates" << flush;

                    //_has_coordinates = true;
                    bool _has_QMAtoms = _orbitals->hasQMAtoms();

                    // three garbage lines
                    getline(_input_file, _line);
                    // now starts the data in format
                    // _id type Qnuc x y z
                    std::vector<std::string> _row;
                    getline(_input_file, _line);
                    boost::trim(_line);

                    boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int nfields = _row.size();

                    int atom_id = 0;
                    while (nfields == 4) {
                        //int atom_id = boost::lexical_cast< int >( _row.at(0) );
                        //int atom_number = boost::lexical_cast< int >( _row.at(0) );
                        std::string _atom_type = _row.at(0);
                        double _x = boost::lexical_cast<double>(_row.at(1));
                        double _y = boost::lexical_cast<double>(_row.at(2));
                        double _z = boost::lexical_cast<double>(_row.at(3));
                        //if ( tools::globals::verbose ) cout << "... ... " << atom_id << " " << atom_type << " " << atom_charge << endl;
                        getline(_input_file, _line);
                        boost::trim(_line);
                        boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                        nfields = _row.size();
                        
                        tools::vec pos=tools::vec(_x,_y,_z);
                        pos*=tools::conv::ang2bohr;

                        if (_has_QMAtoms == false) {
                            _orbitals->AddAtom(atom_id,_atom_type, pos);
                        } else {
                            QMAtom* pAtom = _orbitals->QMAtoms().at(atom_id);
                            pAtom->setPos(pos);
                        }
                        atom_id++;
                    }

                }

                std::string::size_type energy_pos = _line.find("Total Energy");
                if (energy_pos != std::string::npos) {
                    //cout << _line << endl;
                    boost::algorithm::split(results, _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
                    std::string _energy = results[3];
                    boost::trim(_energy);
                    //cout << _energy << endl;
                    _orbitals->setQMEnergy(_conv_Hrt_eV * boost::lexical_cast<double>(_energy));
                    CTP_LOG(ctp::logDEBUG, *_pLog) << (boost::format("QM energy[eV]: %4.6f ") % _orbitals->getQMEnergy()).str() << flush;
                    // _orbitals->_has_qm_energy = true;
                }

                /* Check for ScaHFX = factor of HF exchange included in functional */
                std::string::size_type HFX_pos = _line.find("Fraction HF Exchange ScalHFX");
                if (HFX_pos != std::string::npos) {
                    boost::algorithm::split(results, _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
                    double _ScaHFX = boost::lexical_cast<double>(results.back());
                    _orbitals->setScaHFX(_ScaHFX);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "DFT with " << _ScaHFX << " of HF exchange!" << flush;
                }

                //Finding Basis Dimension, the number of energy levels
                std::string::size_type dim_pos = _line.find("Basis Dimension");
                if (dim_pos != std::string::npos) {

                    boost::algorithm::split(results, _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
                    //_has_basis_dim = true;
                    std::string _dim = results[4]; //The 4th element of results vector is the Basis Dim
                    boost::trim(_dim);
                    _levels = boost::lexical_cast<int>(_dim);
                    //cout <<  boost::lexical_cast<int>(_dim) << endl;
                    //_basis_size = _levels;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Basis Dimension: " << _levels << flush;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Energy levels: " << _levels << flush;
                }
                /********************************************************/


                std::string::size_type OE_pos = _line.find("ORBITAL ENERGIES");
                if (OE_pos != std::string::npos) {
                    getline(_input_file, _line);
                    getline(_input_file, _line);
                    getline(_input_file, _line);
                    if (_line.find("E(Eh)") == std::string::npos) {
                        CTP_LOG(ctp::logDEBUG, *_pLog) << "Warning: Orbital Energies not found in log file" << flush;
                    }
                    for (unsigned i = 0; i < _levels; i++) {
                        getline(_input_file, _line);
                        boost::trim(_line);
                        boost::algorithm::split(results, _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);


                        std::string _no = results[0];

                        boost::trim(_no);
                        unsigned levelnumber = boost::lexical_cast<unsigned>(_no);
                        if (levelnumber != i) {
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "Have a look at the orbital energies something weird is going on" << flush;
                        }
                        std::string _oc = results[1];
                        boost::trim(_oc);
                        double occ = boost::lexical_cast<double>(_oc);
                        // We only count alpha electrons, each orbital must be empty or doubly occupied
                        if (occ == 2 || occ == 1) {
                            _number_of_electrons++;
                            _occ[i] = occ;
                        } else if (occ == 0) {
                            _occ[i] = occ;
                        } else {
                            if (occ == 1){
                                CTP_LOG(ctp::logDEBUG, *_pLog) << "Watch out! No distinction between alpha and beta electrons. Check if occ = 1 is suitable for your calculation " << flush;
                                _number_of_electrons++;
                                _occ[i] = occ;
                            } else {
                            throw runtime_error("Only empty or doubly occupied orbitals are allowed not running the right kind of DFT calculation");
                            }
                        }

                        std::string _e = results[2];
                        boost::trim(_e);
                        _energies [i] = boost::lexical_cast<double>(_e);
                    }
                }
                /*
                 *  Partial charges from the input file
                 */
                std::string::size_type charge_pos = _line.find("CHELPG Charges");

                if (charge_pos != std::string::npos && _get_charges) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting charges" << flush;
                    getline(_input_file, _line);

                    std::vector<std::string> _row;
                    getline(_input_file, _line);
                    boost::trim(_line);
                    boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int nfields = _row.size();

                    while (nfields == 4) {
                        int atom_id = boost::lexical_cast< int >(_row.at(0));
                        atom_id++;
                        std::string atom_type = _row.at(1);
                        double atom_charge = boost::lexical_cast< double >(_row.at(3));
                        getline(_input_file, _line);
                        boost::trim(_line);
                        boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                        nfields = _row.size();
                        
                        QMAtom* pAtom;
                        if (!_orbitals->hasQMAtoms()) {
                            pAtom =_orbitals->AddAtom(atom_id - 1,atom_type, 0, 0, 0);
                        } else {
                            pAtom = _orbitals->QMAtoms().at(atom_id - 1);
                        }
                        pAtom->setPartialcharge(atom_charge);
                    }
                }
            }


            CTP_LOG(ctp::logDEBUG, *_pLog) << "Alpha electrons: " << _number_of_electrons << flush;
            int _occupied_levels = _number_of_electrons;
            int _unoccupied_levels = _levels - _occupied_levels;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Occupied levels: " << _occupied_levels << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Unoccupied levels: " << _unoccupied_levels << flush;


            /************************************************************/

            // copying information to the orbitals object
            /* _orbitals->setBasisSetSize(  _basis_size );*/
            _orbitals->setBasisSetSize(_levels);

            _orbitals->setNumberOfElectrons(_number_of_electrons);

            _orbitals->setNumberOfLevels(_occupied_levels, _unoccupied_levels);

            _orbitals->setSelfEnergy(0.0);

            // copying energies to a matrix
            _orbitals->MOEnergies().resize(_levels);
            //_level = 1;
            for (int i = 0; i < _orbitals->MOEnergies().size(); i++) {
                _orbitals->MOEnergies()[i] = _energies[ i ];
            }


            // cleanup
            // _coefficients.clear();
            _energies.clear();
            _occ.clear();

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done reading Log file" << flush;

            return true;
        }//ParseOrbitalFile(Orbital* _orbital)

        bool Orca::CheckLogFile() {

            // check if the log file exists
            ifstream _input_file((_run_dir + "/" + _log_file_name).c_str());

            if (_input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "Orca LOG is not found" << flush;
                return false;
            };

            std::string _line;
            while (_input_file) {
                getline(_input_file, _line);
                boost::trim(_line);


                std::string::size_type error = _line.find("FATAL ERROR ENCOUNTERED");

                if (error != std::string::npos) {
                    CTP_LOG(ctp::logERROR, *_pLog) << "ORCA encountered a fatal error, maybe a look in the log file may help." << flush;
                    return false;
                }
                error = _line.find("mpirun detected that one or more processes exited with non-zero status");

                if (error != std::string::npos) {
                    CTP_LOG(ctp::logERROR, *_pLog) << "ORCA had an mpi problem, maybe your openmpi version is not good." << flush;
                    return false;
                }
            }
            return true;
        }

        // Parses the Orca gbw file and stores data in the Orbitals object

        bool Orca::ParseOrbitalsFile(Orbitals* _orbitals) {
            if (!CheckLogFile()) return false;
            std::vector<double> _coefficients;
            int _basis_size = _orbitals->getBasisSetSize();
            int _levels = _orbitals->getNumberOfLevels();

            if (_basis_size == 0 || _levels == 0) {
                throw runtime_error("Basis size not set, calculator does not parse log file first");
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Reading the gbw file, this may or may not work so be careful: " << flush;
            ifstream infile;
            infile.open((_run_dir + "/" + _orb_file_name).c_str(), ios::binary | ios::in);
            if (!infile) {
                throw runtime_error("Could not open " + _orb_file_name + " file");
            }
            infile.seekg(24, ios::beg);
            char* buffer = new char [8];
            infile.read(buffer, 8);
            long int offset = *((long int*) buffer);

            infile.seekg(offset, ios::beg);
            infile.read(buffer, 4);
            int op_read = *((int*) buffer);
            infile.seekg(offset + 4, ios::beg);
            infile.read(buffer, 4);
            int dim_read = *((int*) buffer);
            infile.seekg(offset + 8, ios::beg);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Number of operators: " << op_read << " Basis dimension: " << dim_read << flush;
            int n = op_read * dim_read*dim_read;
            delete[] buffer;
            buffer = new char [8];
            for (int i = 0; i < n; i++) {
                infile.read(buffer, 8);
                double mocoeff = *((double*) buffer);
                //CTP_LOG(logDEBUG,*_pLog) << mocoeff<< flush ;
                _coefficients.push_back(mocoeff);
            }
            delete[] buffer;

            infile.close();
            //cout<< "basissize " <<_basis_size << endl;
            //cout << "coeffvektor size "<< _coefficients.size() << endl;


            // i -> MO, j -> AO
            (_orbitals->MOCoefficients()).resize(_levels, _basis_size);
            for (int i = 0; i < _orbitals->MOCoefficients().rows(); i++) {
                for (int j = 0; j < _orbitals->MOCoefficients().cols(); j++) {
                    _orbitals->MOCoefficients()(j, i) = _coefficients[j * _basis_size + i];
                   
                }
            }
           
           ReorderOutput(_orbitals);


            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
            return true;
        }

        std::string Orca::getLName(int lnum) {
            if (lnum == 0) {
                return "S";
            } else if (lnum == 1) {
                return "P";
            } else if (lnum == 2) {
                return "D";
            } else if (lnum == 3) {
                return "F";
            } else {
                throw runtime_error("Orca::getLName functions higher than F not implemented");
            }
            return "0";
        }

        std::string Orca::indent(const double &number) {
            std::stringstream _ssnumber;
            if (number >= 0) {
                _ssnumber << "    ";
            } else {
                _ssnumber << "   ";
            }

            _ssnumber << setiosflags(ios::fixed) << setprecision(15) << std::scientific << number;
            std::string _snumber = _ssnumber.str();
            //boost::replace_first(_snumber, "e", "D");
            return _snumber;
        }




    }
}
