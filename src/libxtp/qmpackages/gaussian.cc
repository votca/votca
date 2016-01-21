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

#include "gaussian.h"
#include "votca/xtp/segment.h"

#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include <stdio.h>
#include <iomanip>
#include <sys/stat.h>
#include <vector>

using namespace std;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        void Gaussian::Initialize(Property *options) {

            // GAUSSIAN file names
            string fileName = "system";

            _xyz_file_name = fileName + ".xyz";
            _input_file_name = fileName + ".com";
            _log_file_name = fileName + ".log";
            _shell_file_name = fileName + ".sh";
            _orb_file_name = "fort.7";
            _input_vxc_file_name = fileName + "-2.com";


            string key = "package";
            string _name = options->get(key + ".name").as<string> ();

            if (_name != "gaussian") {
                cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error("Wrong options file");
            }

            _executable = options->get(key + ".executable").as<string> ();
            _charge = options->get(key + ".charge").as<int> ();
            _spin = options->get(key + ".spin").as<int> ();
            _options = options->get(key + ".options").as<string> ();
            _memory = options->get(key + ".memory").as<string> ();
            _threads = options->get(key + ".threads").as<int> ();
            _chk_file_name = options->get(key + ".checkpoint").as<string> ();
            _scratch_dir = options->get(key + ".scratch").as<string> ();
            _cleanup = options->get(key + ".cleanup").as<string> ();
            
            
            if (options->exists(key + ".vdWRadii")) {
                _vdWfooter = options->get(key + ".vdWRadii").as<string> ();   
            }
            else _vdWfooter="";
            
            
             if (options->exists(key + ".outputVxc")) {
                _output_Vxc = options->get(key + "outputVxc").as<bool> ();   
            }
             else _output_Vxc=false;

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
            iop_pos = _options.find("charge");
            if (iop_pos != std::string::npos) {
                _get_self_energy = true;
                _write_charges = true;
            } else {
                _get_self_energy = false;
                _write_charges = false;
            }

            // check if the basis set is available ("/gen")
            iop_pos = _options.find("gen");
            if (iop_pos != std::string::npos) {
                _write_basis_set = true;
                _basisset_name = options->get(key + ".basisset").as<string> ();
            } else {
                _write_basis_set = false;
            }

            // check if pseudopotentials are required ("pseudo")
            iop_pos = _options.find("pseudo");
            if (iop_pos != std::string::npos) {
                _write_pseudopotentials = true;
            } else {
                _write_pseudopotentials = false;
            }

        }

        /**
         * Prepares the com file from a vector of segments
         * Appends a guess constructed from monomer orbitals if supplied
         */
        bool Gaussian::WriteInputFile(vector<Segment* > segments, Orbitals* orbitals_guess) {
            vector< Atom* > _atoms;
            vector< Atom* > ::iterator ait;
            vector< Segment* >::iterator sit;
            string temp_suffix = "/id";
            string scratch_dir_backup = _scratch_dir;
            //int qmatoms = 0;

            ofstream _com_file;

            string _com_file_name_full = _run_dir + "/" + _input_file_name;

            _com_file.open(_com_file_name_full.c_str());
            // header 
            if (_chk_file_name.size()) _com_file << "%chk=" << _chk_file_name << endl;
            if (_memory.size()) _com_file << "%mem=" << _memory << endl;
            if (_threads > 0) _com_file << "%nprocshared=" << _threads << endl;
            if (_options.size()) _com_file << _options << endl;

            _com_file << endl;
            _com_file << "TITLE ";

            for (sit = segments.begin(); sit != segments.end(); ++sit) {
                _com_file << (*sit)->getName() << " ";
            }
            _com_file << endl << endl;
            _com_file << setw(2) << _charge << setw(2) << _spin << endl;

            // write coordinates!
            if (!_write_charges) {

                for (sit = segments.begin(); sit != segments.end(); ++sit) {
                    temp_suffix = temp_suffix + "_" + boost::lexical_cast<string>((*sit)->getId());
                    _atoms = (*sit)-> Atoms();

                    for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                        if ((*ait)->HasQMPart() == false) {
                            continue;
                        }

                        vec pos = (*ait)->getQMPos();
                        string name = (*ait)->getElement();

                        //fprintf(out, "%2s %4.7f %4.7f %4.7f \n"
                        _com_file << setw(3) << name.c_str()
                                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX()*10
                                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY()*10
                                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ()*10
                                << endl;
                    }
                }

                // BASISSET
                if (_write_basis_set) {

                    _com_file << endl;
                    list<string> elements;
                    BasisSet bs;
                    // string basis_name(_basis);

                    bs.LoadBasisSet(_basisset_name);
                    LOG(logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;

                    for (sit = segments.begin(); sit != segments.end(); ++sit) {

                        vector< Atom* > atoms = (*sit)-> Atoms();
                        vector< Atom* >::iterator it;

                        for (it = atoms.begin(); it < atoms.end(); it++) {

                            string element_name = (*it)->getElement();

                            list<string>::iterator ite;
                            ite = find(elements.begin(), elements.end(), element_name);

                            if (ite == elements.end()) {
                                elements.push_back(element_name);

                                Element* element = bs.getElement(element_name);
                                /* Alternative is to write each basis set to a element_name.gbs file
                                 * and include the gbs file in the com-file via Gaussian's @ function
                                 * Advantage: *gbs files can be reused by isogwa later
                                 */
                                ofstream _el_file;
                                string _el_file_name = _run_dir + "/" + element_name + ".gbs";
                                _el_file.open(_el_file_name.c_str());
                                // element name, [possibly indeces of centers], zero to indicate the end
                                //_com_file << element_name << " 0" << endl;
                                _com_file << "@" << element_name << ".gbs" << endl;
                                _el_file << element_name << " 0" << endl;
                                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {

                                    Shell* shell = (*its);
                                    // shell type, number primitives, scale factor
                                    //_com_file << shell->getType() << " " << shell->getSize() << " " << shell->getScale() << endl;
                                    _el_file << shell->getType() << " " << shell->getSize() << " " << FortranFormat(shell->getScale()) << endl;
                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                        //_com_file << gaussian->decay << " " << gaussian->contraction << endl;
                                        _el_file << FortranFormat(gaussian->decay);
                                        for (unsigned _icontr = 0; _icontr < gaussian->contraction.size(); _icontr++) {
                                            if (gaussian->contraction[_icontr] != 0.0) {
                                                _el_file << " " << FortranFormat(gaussian->contraction[_icontr]);
                                            }
                                        }
                                        _el_file << endl;
                                    }
                                }

                                //_com_file << "****\n";
                                _el_file << "****\n";
                                _el_file.close();

                            }
                        }
                    }
                } // end BASISSET


                // ECP
                if (_write_pseudopotentials) {
                    string pseudopotential_name("ecp");

                    _com_file << endl;
                    list<string> elements;

                    elements.push_back("H");
                    elements.push_back("He");

                    BasisSet ecp;
                    ecp.LoadPseudopotentialSet(pseudopotential_name);

                    LOG(logDEBUG, *_pLog) << "Loaded Pseudopotentials " << pseudopotential_name << flush;

                    for (sit = segments.begin(); sit != segments.end(); ++sit) {

                        vector< Atom* > atoms = (*sit)-> Atoms();
                        vector< Atom* >::iterator it;

                        for (it = atoms.begin(); it < atoms.end(); it++) {

                            string element_name = (*it)->getElement();

                            list<string>::iterator ite;
                            ite = find(elements.begin(), elements.end(), element_name);

                            if (ite == elements.end()) {
                                elements.push_back(element_name);

                                Element* element = ecp.getElement(element_name);
                                   
                                // element name, [possibly indeces of centers], zero to indicate the end
                                _com_file << element_name << " 0\n"
                                        << pseudopotential_name << " "
                                        << element->getLmax() << " " << element->getNcore() << endl;
                                //write local component
                                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {

                                    Shell* shell = (*its);
                                    // shell type, number primitives, scale factor
                                    if (shell->getLmax() == element->getLmax()) {
                                        _com_file << shell->getType() << endl;
                                        _com_file << shell->getSize() << endl;

                                        for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                            GaussianPrimitive* gaussian = *itg;
                                            _com_file << gaussian->power << " " << gaussian->decay << " " << gaussian->contraction[0] << endl;
                                        }
                                    }
                                }
                                // write remaining shells in ascending order s,p,d...
                                for (int i = 0; i < element->getLmax(); i++) {
                                    for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                        Shell* shell = (*its);
                                        if (shell->getLmax() == i) {
                                            // shell type, number primitives, scale factor
                                            _com_file << shell->getType() << endl;
                                            _com_file << shell->getSize() << endl;

                                            for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                                GaussianPrimitive* gaussian = *itg;
                                                _com_file << gaussian->power << " " << gaussian->decay << " " << gaussian->contraction[0] << endl;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }// end ECP

                





            } else {

                vector< QMAtom* > *qmatoms = orbitals_guess->getAtoms();
                vector< QMAtom* >::iterator it;

                // This is needed for the QM/MM scheme, since only orbitals have 
                // updated positions of the QM region, hence vector<Segments*> is 
                // NULL in the QMMachine and the QM region is also printed here
                for (it = qmatoms->begin(); it < qmatoms->end(); it++) {
                    if (!(*it)->from_environment) {
                        _com_file << setw(3) << (*it)->type.c_str()
                                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << (*it)->x
                                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << (*it)->y
                                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << (*it)->z
                                << endl;


                        //_com_file << (*it)->type << " " <<  (*it)->x << " " << (*it)->y << " " << (*it)->z << endl;
                    }
                }

                _com_file << endl;

                // if we need to write basis sets, do it now
                if (_write_basis_set) {

                    list<string> elements;
                    BasisSet bs;
                    // string basis_name(_basis);

                    bs.LoadBasisSet(_basisset_name);
                    LOG(logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;

                    for (it = qmatoms->begin(); it < qmatoms->end(); it++) {
                        if (!(*it)->from_environment) {
                            string element_name = (*it)->type;

                            //cout << "looking up basis set for element " << element_name << endl;

                            list<string>::iterator ite;
                            ite = find(elements.begin(), elements.end(), element_name);

                            if (ite == elements.end()) {
                                elements.push_back(element_name);

                                Element* element = bs.getElement(element_name);
                                /* Alternative is to write each basis set to a element_name.gbs file
                                 * and include the gbs file in the com-file via Gaussian's @ function
                                 * Advantage: *gbs files can be reused by isogwa later
                                 */
                                ofstream _el_file;
                                string _el_file_name = _run_dir + "/" + element_name + ".gbs";
                                _el_file.open(_el_file_name.c_str());
                                // element name, [possibly indeces of centers], zero to indicate the end
                                //_com_file << element_name << " 0" << endl;
                                _com_file << "@" << element_name << ".gbs" << endl;
                                _el_file << element_name << " 0" << endl;
                                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {

                                    Shell* shell = (*its);
                                    // shell type, number primitives, scale factor
                                    //_com_file << shell->getType() << " " << shell->getSize() << " " << shell->getScale() << endl;
                                    _el_file << shell->getType() << " " << shell->getSize() << " " << FortranFormat(shell->getScale()) << endl;
                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                        //_com_file << gaussian->decay << " " << gaussian->contraction << endl;
                                        _el_file << FortranFormat(gaussian->decay);
                                        for (unsigned _icontr = 0; _icontr < gaussian->contraction.size(); _icontr++) {
                                            if (gaussian->contraction[_icontr] != 0.0) {
                                                _el_file << " " << FortranFormat(gaussian->contraction[_icontr]);
                                            }
                                        }
                                        _el_file << endl;
                                    }
                                }

                                //_com_file << "****\n";
                                _el_file << "****\n";
                                _el_file.close();

                            }
                        }
                    }

                }

                if (_write_pseudopotentials) {
                    string pseudopotential_name("ecp");

                    _com_file << endl;
                    list<string> elements;

                    elements.push_back("H");
                    elements.push_back("He");

                    BasisSet ecp;
                    ecp.LoadPseudopotentialSet(pseudopotential_name);

                    LOG(logDEBUG, *_pLog) << "Loaded Pseudopotentials " << pseudopotential_name << flush;

                    //for (sit = segments.begin(); sit != segments.end(); ++sit) {

                    //  vector< Atom* > atoms = (*sit)-> Atoms();
                    // vector< Atom* >::iterator it;

                    for (it = qmatoms->begin(); it < qmatoms->end(); it++) {
                        if (!(*it)->from_environment) {
                            string element_name = (*it)->type;

                            list<string>::iterator ite;
                            ite = find(elements.begin(), elements.end(), element_name);

                            if (ite == elements.end()) {
                                elements.push_back(element_name);

                                Element* element = ecp.getElement(element_name);

                                // element name, [possibly indeces of centers], zero to indicate the end
                                _com_file << element_name << " 0\n"
                                        << pseudopotential_name << " "
                                        << element->getLmax() << " " << element->getNcore() << endl;

                                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {

                                    Shell* shell = (*its);
                                    // shell type, number primitives, scale factor
                                    _com_file << shell->getType() << endl;
                                    _com_file << shell->getSize() << endl;

                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                        _com_file << gaussian->power << " " << FortranFormat(gaussian->decay) << " " << FortranFormat(gaussian->contraction[0]) << endl;
                                    }
                                }
                            }
                        }
                    }
                    // }
                    _com_file << endl;
                } // ECP


                // actually write the charges

                for (it = qmatoms->begin(); it < qmatoms->end(); it++) {
                    if ((*it)->from_environment) {
                        boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
                        fmt % (*it)->x % (*it)->y % (*it)->z % (*it)->charge;
                        if ((*it)->charge != 0.0) _com_file << fmt << endl;
                    }
                }

                _com_file << endl;
            } // end writing coordinates depending on atoms or qmatoms presence

            if (_write_guess) {
                if (orbitals_guess == NULL) {
                    throw std::runtime_error("A guess for dimer orbitals has not been prepared.");
                } else {
                    vector<int> _sort_index;

                    orbitals_guess->SortEnergies(&_sort_index);

                    _com_file << endl << "(5D15.8)" << endl;

                    int level = 1;
                    int ncolumns = 5;

                    for (vector< int > ::iterator soi = _sort_index.begin(); soi != _sort_index.end(); ++soi) {

                        double _energy = (orbitals_guess->_mo_energies)[*soi];

                        _com_file << setw(5) << level << " Alpha MO OE=" << FortranFormat(_energy) << endl;

                        ub::matrix_row< ub::matrix<double> > mr(orbitals_guess->_mo_coefficients, *soi);

                        int column = 1;
                        for (unsigned j = 0; j < mr.size(); ++j) {
                            _com_file << FortranFormat(mr[j]);
                            if (column == ncolumns) {
                                _com_file << std::endl;
                                column = 0;
                            }
                            column++;
                        }

                        level++;
                        if (column != 1) _com_file << endl;
                    }
                }
            }





            if (_output_Vxc) {
             

                ofstream _com_file2;

                string _com_file_name_full2 = _run_dir + "/" + _input_vxc_file_name;

                _com_file2.open(_com_file_name_full2.c_str());
                // header 
                if (_chk_file_name.size()) _com_file2 << "%chk=" << _chk_file_name << endl;
                if (_memory.size()) _com_file2 << "%mem=" << _memory << endl;
                _com_file2 << "%nprocshared=1" << endl;

                // adjusting the options line to Vxc output only
                string _options_vxc = _options;
                boost::algorithm::replace_all(_options_vxc, "pseudo=read", "Geom=AllCheck");
                boost::algorithm::replace_all(_options_vxc, "/gen", " chkbasis");
                boost::algorithm::replace_all(_options_vxc, "punch=mo", "guess=read");
                boost::algorithm::replace_all(_options_vxc, "guess=tcheck", "");
                boost::algorithm::replace_all(_options_vxc, "guess=huckel", "");
                boost::algorithm::replace_all(_options_vxc, "charge", "charge=check");
                if (_options_vxc.size()) _com_file2 << _options_vxc << endl;

                // # pop=minimal pbepbe/gen pseudo=read scf=tight punch=mo
                // # pop=minimal pbepbe chkbasis nosymm Geom=AllCheck Guess=Read


                _com_file2 << endl;
                _com_file2 << "VXC output run \n";
                _com_file2 << endl;
                _com_file2.close();




            }

            _com_file << _vdWfooter << endl;


            _com_file << endl;
            _com_file.close();
            // and now generate a shell script to run both jobs
            LOG(logDEBUG, *_pLog) << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

            _scratch_dir = scratch_dir_backup + temp_suffix;
            
            //boost::filesystem::create_directories(_scratch_dir + temp_suffix);
            //string _temp("scratch_dir " + _scratch_dir + temp_suffix + "\n");
            //_com_file << _temp;
            WriteShellScript();
            _scratch_dir = scratch_dir_backup;
            
            
            return true;
        }

        bool Gaussian::WriteShellScript() {
            ofstream _shell_file;

            string _shell_file_name_full = _run_dir + "/" + _shell_file_name;

            _shell_file.open(_shell_file_name_full.c_str());

            _shell_file << "#!/bin/tcsh" << endl;
            _shell_file << "mkdir -p " << _scratch_dir << endl;
            _shell_file << "setenv GAUSS_SCRDIR " << _scratch_dir << endl;
            _shell_file << _executable << " " << _input_file_name << endl;
            if (_output_Vxc) {
                _shell_file << "rm fort.22" << endl;
                _shell_file << "setenv DoPrtXC YES" << endl;
                _shell_file << _executable << " " << _input_vxc_file_name << " >& /dev/null " << endl;
                _shell_file << "setenv DoPrtXC NO" << endl;
                _shell_file << "rm $GAUSS_SCRDIR/*" << endl;
            }
            _shell_file.close();

            return true;
        }

        /**
         * Runs the Gaussian job. 
         */
        bool Gaussian::Run() {

            LOG(logDEBUG, *_pLog) << "GAUSSIAN: running [" << _executable << " " << _input_file_name << "]" << flush;

            if (system(NULL)) {
                // if scratch is provided, run the shell script; 
                // otherwise run gaussian directly and rely on global variables 
                string _command;
                if (_scratch_dir.size() != 0 || _output_Vxc) {
                    _command = "cd " + _run_dir + "; tcsh " + _shell_file_name;
                    //            _command  = "cd " + _run_dir + "; mkdir -p " + _scratch_dir +"; " + _executable + " " + _input_file_name;
                } else {
                    _command = "cd " + _run_dir + "; mkdir -p $GAUSS_SCRDIR; " + _executable + " " + _input_file_name;
                }

                //int i = system(_command.c_str());
                system(_command.c_str());
                if (CheckLogFile()) {
                    LOG(logDEBUG, *_pLog) << "GAUSSIAN: finished job" << flush;
                    return true;
                } else {
                    LOG(logDEBUG, *_pLog) << "GAUSSIAN: job failed" << flush;
                }
            } else {
                LOG(logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
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

                LOG(logDEBUG, *_pLog) << "Removing " << _cleanup << " files" << flush;
                Tokenizer tok_cleanup(_cleanup, ",");
                vector <string> _cleanup_info;
                tok_cleanup.ToVector(_cleanup_info);

                vector<string> ::iterator it;

                for (it = _cleanup_info.begin(); it != _cleanup_info.end(); ++it) {

                    if (*it == "com") {
                        string file_name = _run_dir + "/" + _input_file_name;
                        remove(file_name.c_str());
                    }

                    if (*it == "sh") {
                        string file_name = _run_dir + "/" + _shell_file_name;
                        remove(file_name.c_str());
                    }

                    if (*it == "log") {
                        string file_name = _run_dir + "/" + _log_file_name;
                        remove(file_name.c_str());
                    }

                    if (*it == "chk") {
                        string file_name = _run_dir + "/" + _chk_file_name;
                        remove(file_name.c_str());
                    }

                    if (*it == "fort.7") {
                        string file_name = _run_dir + "/" + *it;
                        remove(file_name.c_str());
                    }
                }
            }

        }

        /**
         * Reads in the MO coefficients from a GAUSSIAN fort.7 file
         */
        bool Gaussian::ParseOrbitalsFile(Orbitals * _orbitals) {
            std::map <int, std::vector<double> > _coefficients;
            std::map <int, double> _energies;

            std::string _line;
            unsigned _levels = 0;
            unsigned _level;
            unsigned _basis_size = 0;

            string _orb_file_name_full = _orb_file_name;
            if (_run_dir != "") _orb_file_name_full = _run_dir + "/" + _orb_file_name;
            std::ifstream _input_file(_orb_file_name_full.c_str());

            if (_input_file.fail()) {
                LOG(logERROR, *_pLog) << "File " << _orb_file_name << " with molecular orbitals is not found " << flush;
                return false;
            } else {
                LOG(logDEBUG, *_pLog) << "Reading MOs from " << _orb_file_name << flush;
            }

            // number of coefficients per line is  in the first line of the file (5D15.8)
            getline(_input_file, _line);
            std::vector<string> strs;
            boost::algorithm::split(strs, _line, boost::is_any_of("(D)"));
            //clog << strs.at(1) << endl;
            //int nrecords_in_line = boost::lexical_cast<int>(strs.at(1));
            string format = strs.at(2);

            //clog << endl << "Orbital file " << filename << " has " 
            //        << nrecords_in_line << " records per line, in D"
            //        << format << " format." << endl;

            while (_input_file) {

                getline(_input_file, _line);
                // if a line has an equality sign, must be energy
                std::string::size_type energy_pos = _line.find("=");

                if (energy_pos != std::string::npos) {

                    vector<string> results;
                    boost::trim(_line);

                    boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                            boost::algorithm::token_compress_on);
                    //cout << results[1] << ":" << results[2] << ":" << results[3] << ":" << results[4] << endl;

                    _level = boost::lexical_cast<int>(results.front());
                    boost::replace_first(results.back(), "D", "e");
                    _energies[ _level ] = boost::lexical_cast<double>(results.back());
                    _levels++;

                } else {

                    while (_line.size() > 1) {
                        string _coefficient;
                        _coefficient.assign(_line, 0, 15);
                        boost::trim(_coefficient);
                        boost::replace_first(_coefficient, "D", "e");
                        double coefficient = boost::lexical_cast<double>(_coefficient);
                        _coefficients[ _level ].push_back(coefficient);
                        _line.erase(0, 15);
                    }
                }
            }

            // some sanity checks
            LOG(logDEBUG, *_pLog) << "Energy levels: " << _levels << flush;

            std::map< int, vector<double> >::iterator iter = _coefficients.begin();
            _basis_size = iter->second.size();

            for (iter = _coefficients.begin()++; iter != _coefficients.end(); iter++) {
                if (iter->second.size() != _basis_size) {
                    LOG(logERROR, *_pLog) << "Error reading " << _orb_file_name << ". Basis set size change from level to level." << flush;
                    return false;
                }
            }

            LOG(logDEBUG, *_pLog) << "Basis set size: " << _basis_size << flush;

            // copying information to the orbitals object
            _orbitals->setBasisSetSize(_basis_size); // = _basis_size;
            // _orbitals->_has_basis_set_size = true;
            // _orbitals->_has_mo_coefficients = true;
            // _orbitals->_has_mo_energies = true;

            // copying energies to the orbitals object  
            ub::vector<double> &mo_energies = _orbitals->MOEnergies();
            mo_energies.resize(_levels);
            for (size_t i = 0; i < mo_energies.size(); i++) mo_energies[i] = _energies[ i + 1 ];

            // copying mo coefficients to the orbitals object
            ub::matrix<double> &mo_coefficients = _orbitals->MOCoefficients();
            mo_coefficients.resize(_levels, _basis_size);
            for (size_t i = 0; i < mo_coefficients.size1(); i++)
                for (size_t j = 0; j < mo_coefficients.size2(); j++)
                    _orbitals->_mo_coefficients(i, j) = _coefficients[i + 1][j];


            //cout << _mo_energies << endl;   
            //cout << _mo_coefficients << endl; 

            LOG(logDEBUG, *_pLog) << "GAUSSIAN: done reading MOs" << flush;

            return true;
        }

        bool Gaussian::CheckLogFile() {

            // check if the log file exists
            boost::filesystem::path arg_path;
            char ch;

            string _full_name = (arg_path / _run_dir / _log_file_name).c_str();
            ifstream _input_file(_full_name.c_str());

            if (_input_file.fail()) {
                LOG(logERROR, *_pLog) << "GAUSSIAN: " << _full_name << " is not found" << flush;
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
            } while (ch != '\n' && (int) _input_file.tellg() != -1);

            string _line;
            getline(_input_file, _line); // Read the current line
            //cout << "\nResult: " << _line << '\n';     // Display it
            _input_file.close();

            std::string::size_type self_energy_pos = _line.find("Normal termination of Gaussian");
            if (self_energy_pos == std::string::npos) {
                LOG(logERROR, *_pLog) << "GAUSSIAN: " << _full_name << " is incomplete" << flush;
                return false;
            } else {
                //LOG(logDEBUG,*_pLog) << "Gaussian LOG is complete" << flush;
                return true;
            }
        }

        /**
         * Parses the Gaussian Log file and stores data in the Orbitals object 
         */
        bool Gaussian::ParseLogFile(Orbitals * _orbitals) {

            static const double _conv_Hrt_eV = 27.21138386;

            string _line;
            vector<string> results;
            bool _has_occupied_levels = false;
            bool _has_unoccupied_levels = false;
            bool _has_number_of_electrons = false;
            bool _has_basis_set_size = false;
            bool _has_overlap_matrix = false;
            //bool _has_vxc_matrix = false;
            bool _has_charges = false;
            //bool _has_coordinates = false;
            //bool _has_qm_energy = false;
            bool _has_self_energy = false;

            bool _read_vxc = false;

            int _occupied_levels = 0;
            int _unoccupied_levels = 0;
            int _number_of_electrons = 0;
            int _basis_set_size = 0;
            int _cart_basis_set_size = 0;

            LOG(logDEBUG, *_pLog) << "GAUSSIAN: parsing " << _log_file_name << flush;

            string _log_file_name_full = _log_file_name;
            if (_run_dir != "") _log_file_name_full = _run_dir + "/" + _log_file_name;

            // check if LOG file is complete
            if (!CheckLogFile()) return false;

            // save qmpackage name
            //_orbitals->_has_qm_package = true;
            _orbitals->setQMpackage("gaussian");


            // Start parsing the file line by line
            ifstream _input_file(_log_file_name_full.c_str());
            while (_input_file) {

                getline(_input_file, _line);
                boost::trim(_line);

                /*
                 * Check is pseudo keyword is present in LOG file -> read vxc
                 */
                std::string::size_type pseudo_pos = _line.find("pseudo=read");
                if (pseudo_pos != std::string::npos) {
                    _read_vxc = true;
                    // Uncomment for next version
                    //_orbitals->setWithECP(true);
                }

                /* Check for ScaHFX = factor of HF exchange included in functional */
                std::string::size_type HFX_pos = _line.find("ScaHFX=");
                if (HFX_pos != std::string::npos) {
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    double _ScaHFX = boost::lexical_cast<double>(results.back());
                    _orbitals->setScaHFX(_ScaHFX);
                    LOG(logDEBUG, *_pLog) << "DFT with " << _ScaHFX << " of HF exchange!" << flush;
                }



                /*
                 * number of occupied and virtual orbitals
                 * N alpha electrons      M beta electrons
                 */
                std::string::size_type electrons_pos = _line.find("alpha electrons");
                if (electrons_pos != std::string::npos) {
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    _has_number_of_electrons = true;
                    _number_of_electrons = boost::lexical_cast<int>(results.front());
                    _orbitals->setNumberOfElectrons(_number_of_electrons);
                    // _orbitals->_has_number_of_electrons = true;
                    LOG(logDEBUG, *_pLog) << "Alpha electrons: " << _number_of_electrons << flush;
                }

                /*
                 * basis set size
                 * N basis functions,  M primitive gaussians,   K cartesian basis functions
                 */
                std::string::size_type basis_pos = _line.find("basis functions,");
                if (basis_pos != std::string::npos) {
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    _has_basis_set_size = true;
                    _basis_set_size = boost::lexical_cast<int>(results.front());
                    _orbitals->setBasisSetSize(_basis_set_size);
                    // _orbitals->_has_basis_set_size = true;
                    _cart_basis_set_size = boost::lexical_cast<int>(results[6]);
                    LOG(logDEBUG, *_pLog) << "Basis functions: " << _basis_set_size << flush;
                    if (_read_vxc) {
                        LOG(logDEBUG, *_pLog) << "Cartesian functions: " << _cart_basis_set_size << flush;
                    }
                }

                /*
                 * energies of occupied/unoccupied levels
                 * Alpha  occ.(virt.) eigenvalues -- e1 e2 e3 e4 e5
                 */
                std::string::size_type eigenvalues_pos = _line.find("Alpha");
                if (eigenvalues_pos != std::string::npos) {

                    std::list<std::string> stringList;
                    //int _unoccupied_levels = 0;
                    //int _occupied_levels = 0;

                    while (eigenvalues_pos != std::string::npos && !_has_occupied_levels && !_has_unoccupied_levels) {
                        //cout << _line << endl;

                        boost::iter_split(stringList, _line, boost::first_finder("--"));

                        vector<string> energies;
                        boost::trim(stringList.back());

                        boost::algorithm::split(energies, stringList.back(), boost::is_any_of("\t "), boost::algorithm::token_compress_on);

                        if (stringList.front().find("virt.") != std::string::npos) {
                            _unoccupied_levels += energies.size();
                            energies.clear();
                        }

                        if (stringList.front().find("occ.") != std::string::npos) {
                            _occupied_levels += energies.size();
                            energies.clear();
                        }

                        getline(_input_file, _line);
                        eigenvalues_pos = _line.find("Alpha");
                        boost::trim(_line);

                        //boost::iter_split(stringList, _line, boost::first_finder("--"));

                        if (eigenvalues_pos == std::string::npos) {
                            _has_occupied_levels = true;
                            _has_unoccupied_levels = true;
                            _orbitals->setNumberOfLevels(_occupied_levels, _unoccupied_levels);
                            // _orbitals->_occupied_levels = _occupied_levels;
                            // _orbitals->_unoccupied_levels = _unoccupied_levels;
                            // _orbitals->_has_occupied_levels = true;
                            // _orbitals->_has_unoccupied_levels = true;
                            LOG(logDEBUG, *_pLog) << "Occupied levels: " << _occupied_levels << flush;
                            LOG(logDEBUG, *_pLog) << "Unoccupied levels: " << _unoccupied_levels << flush;
                        }
                    } // end of the while loop              
                } // end of the eigenvalue parsing


                /*
                 * overlap matrix
                 * stored after the *** Overlap *** line
                 */
                std::string::size_type overlap_pos = _line.find("*** Overlap ***");
                if (overlap_pos != std::string::npos) {

                    // prepare the container
                    ub::symmetric_matrix<double> &overlap = _orbitals->AOOverlap();

                    // _orbitals->_has_overlap = true;
                    overlap.resize(_basis_set_size);

                    _has_overlap_matrix = true;
                    //cout << "Found the overlap matrix!" << endl;   
                    vector<int> _j_indeces;

                    int _n_blocks = 1 + ((_basis_set_size - 1) / 5);
                    //cout << _n_blocks;

                    getline(_input_file, _line);
                    boost::trim(_line);

                    for (int _block = 0; _block < _n_blocks; _block++) {

                        // first line gives the j index in the matrix
                        //cout << _line << endl;

                        boost::tokenizer<> tok(_line);
                        std::transform(tok.begin(), tok.end(), std::back_inserter(_j_indeces), &boost::lexical_cast<int, std::string>);
                        //std::copy( _j_indeces.begin(), _j_indeces.end(), std::ostream_iterator<int>(std::cout,"\n") );

                        // read the block of max _basis_size lines + the following header
                        for (int i = 0; i <= _basis_set_size; i++) {
                            getline(_input_file, _line);
                            //cout << _line << endl;
                            if (std::string::npos == _line.find("D")) break;

                            // split the line on the i index and the rest

                            vector<string> _row;
                            boost::trim(_line);
                            boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);


                            int _i_index = boost::lexical_cast<int>(_row.front());
                            _row.erase(_row.begin());

                            //cout << _i_index << ":" << _line << endl ;

                            std::vector<int>::iterator _j_iter = _j_indeces.begin();

                            for (std::vector<string>::iterator iter = _row.begin()++; iter != _row.end(); iter++) {
                                string _coefficient = *iter;

                                boost::replace_first(_coefficient, "D", "e");
                                //cout << boost::lexical_cast<double>( _coefficient ) << endl;

                                int _j_index = *_j_iter;
                                //_overlap( _i_index-1 , _j_index-1 ) = boost::lexical_cast<double>( _coefficient );
                                overlap(_i_index - 1, _j_index - 1) = boost::lexical_cast<double>(_coefficient);
                                _j_iter++;

                            }


                        }

                        // clear the index for the next block
                        _j_indeces.clear();
                    } // end of the blocks
                    LOG(logDEBUG, *_pLog) << "Read the overlap matrix" << flush;
                } // end of the if "Overlap" found   


                /*
                 *  Partial charges from the input file
                 */
                std::string::size_type charge_pos = _line.find("Charges from ESP fit, RMS");

                if (charge_pos != std::string::npos && _get_charges) {
                    LOG(logDEBUG, *_pLog) << "Getting charges" << flush;
                    _has_charges = true;
                    getline(_input_file, _line);
                    getline(_input_file, _line);

                    bool _has_atoms = _orbitals->hasQMAtoms();

                    vector<string> _row;
                    getline(_input_file, _line);
                    boost::trim(_line);
                    //cout << _line << endl;
                    boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int nfields = _row.size();
                    //cout << _row.size() << endl;

                    while (nfields == 3) {
                        int atom_id = boost::lexical_cast< int >(_row.at(0));
                        //int atom_number = boost::lexical_cast< int >(_row.at(0));
                        string atom_type = _row.at(1);
                        double atom_charge = boost::lexical_cast< double >(_row.at(2));
                        //if ( tools::globals::verbose ) cout << "... ... " << atom_id << " " << atom_type << " " << atom_charge << endl;
                        getline(_input_file, _line);
                        boost::trim(_line);
                        boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                        nfields = _row.size();

                        if (_has_atoms == false) {
                            _orbitals->AddAtom(atom_type, 0, 0, 0, atom_charge);
                        } else {
                            QMAtom* pAtom = _orbitals->_atoms.at(atom_id - 1);
                            pAtom->type = atom_type;
                            pAtom->charge = atom_charge;
                        }

                    }
                    //_orbitals->_has_atoms = true;
                }


                /*
                 * Coordinates of the final configuration
                 * stored in the archive at the end of the file
                 */
                int cpn = 0; // marker appearence marker
                std::string::size_type coordinates_pos = _line.find("\\");

                if (coordinates_pos != std::string::npos && cpn == 0) {
                    ++cpn; // updates but ignores
                    LOG(logDEBUG, *_pLog) << "Getting the coordinates" << flush;
                    //_has_coordinates = true;
                    boost::trim(_line);
                    string archive = _line;
                    while (_line.size() != 0) {
                        getline(_input_file, _line);
                        boost::trim(_line);
                        archive += _line;
                    }

                    bool _has_atoms = _orbitals->hasQMAtoms();
                    std::list<std::string> stringList;
                    vector<string> results;
                    boost::iter_split(stringList, archive, boost::first_finder("\\\\"));

                    list<string>::iterator coord_block = stringList.begin();
                    advance(coord_block, 3);

                    vector<string> atom_block;
                    boost::algorithm::split(atom_block, *coord_block, boost::is_any_of("\\"), boost::algorithm::token_compress_on);

                    vector<string>::iterator atom_block_it;
                    int aindex = 0;

                    for (atom_block_it = ++atom_block.begin(); atom_block_it != atom_block.end(); ++atom_block_it) {
                        vector<string> atom;

                        boost::algorithm::split(atom, *atom_block_it, boost::is_any_of(","), boost::algorithm::token_compress_on);
                        string _atom_type = atom.front();

                        vector<string>::iterator it_atom;
                        it_atom = atom.end();
                        double _z = boost::lexical_cast<double>(*(--it_atom));
                        double _y = boost::lexical_cast<double>(*(--it_atom));
                        double _x = boost::lexical_cast<double>(*(--it_atom));

                        if (_has_atoms == false) {
                            _orbitals->AddAtom(_atom_type, _x, _y, _z);
                        } else {
                            QMAtom* pAtom = _orbitals->_atoms.at(aindex);
                            pAtom->type = _atom_type;
                            pAtom->x = _x;
                            pAtom->y = _y;
                            pAtom->z = _z;
                            aindex++;
                        }

                    }

                    // get the QM energy out
                    advance(coord_block, 1);
                    vector<string> block;
                    vector<string> energy;
                    boost::algorithm::split(block, *coord_block, boost::is_any_of("\\"), boost::algorithm::token_compress_on);
                    //boost::algorithm::split(energy, block[1], boost::is_any_of("="), boost::algorithm::token_compress_on);
                    //_orbitals->setQMEnergy( _conv_Hrt_eV * boost::lexical_cast<double> ( energy[1] ) );
                    map<string, string> properties;
                    vector<string>::iterator block_it;
                    for (block_it = block.begin(); block_it != block.end(); ++block_it) {
                        vector<string> property;
                        boost::algorithm::split(property, *block_it, boost::is_any_of("="), boost::algorithm::token_compress_on);
                        properties[property[0]] = property[1];
                    }
                    LOG(logDEBUG, *_pLog) << "QM energy " << _orbitals->getQMEnergy() << flush;
                    //_has_qm_energy = true;
                    //_orbitals->_has_atoms = true;
                    //_orbitals->_has_qm_energy = true;
                    if (properties.count("HF") > 0) {
                        double energy_hartree = boost::lexical_cast<double>(properties["HF"]);
                        //_orbitals->setQMEnergy(_has_qm_energy = true;
                        _orbitals-> setQMEnergy(_conv_Hrt_eV * energy_hartree);
                        LOG(logDEBUG, *_pLog) << "QM energy " << _orbitals->_qm_energy << flush;
                    } else {
                        cout << endl;
                        throw std::runtime_error("ERROR No energy in archive");
                    }

                    //            boost::algorithm::split(energy, block[1], boost::is_any_of("="), boost::algorithm::token_compress_on);
                    //            cout << endl << energy[1] << endl;
                    //            _orbitals->_qm_energy = _conv_Hrt_eV * boost::lexical_cast<double> ( energy[1] );
                    //            
                    //            LOG(logDEBUG, *_pLog) << "QM energy " << _orbitals->_qm_energy <<  flush;
                    //            _has_qm_energy = true;

                }

                /*
                 * Self-energy of external charges
                 */
                std::string::size_type self_energy_pos = _line.find("Self energy of the charges");

                if (self_energy_pos != std::string::npos) {
                    LOG(logDEBUG, *_pLog) << "Getting the self energy\n";
                    vector<string> block;
                    vector<string> energy;
                    boost::algorithm::split(block, _line, boost::is_any_of("="), boost::algorithm::token_compress_on);
                    boost::algorithm::split(energy, block[1], boost::is_any_of("\t "), boost::algorithm::token_compress_on);

                    // _orbitals->_has_self_energy = true;
                    _orbitals->setSelfEnergy(_conv_Hrt_eV * boost::lexical_cast<double> (energy[1]));

                    LOG(logDEBUG, *_pLog) << "Self energy " << _orbitals->getSelfEnergy() << flush;

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

            LOG(logDEBUG, *_pLog) << "Done parsing" << flush;
            _input_file.close();

            /* Now, again the somewhat ugly construction:
             * if we request writing of pseudopotential data to the input file, this
             * implies a GW-BSE run. For this, we have to 
             * - parse atomic orbitals Vxc matrix */
            if (_read_vxc) {
                LOG(logDEBUG, *_pLog) << "Parsing fort.24 for Vxc" << flush;
                string _log_file_name_full;
                if (_run_dir == "") {
                    _log_file_name_full = "fort.24";
                } else {
                    _log_file_name_full = _run_dir + "/fort.24";
                }


                // prepare the container
                // _orbitals->_has_vxc = true;
                ub::symmetric_matrix<double>& _vxc = _orbitals->AOVxc();
                _vxc.resize(_cart_basis_set_size);


               // _has_vxc_matrix = true;
                //cout << "Found the overlap matrix!" << endl;   
                vector<int> _j_indeces;


                // Start parsing the file line by line
                ifstream _input_file(_log_file_name_full.c_str());
                while (_input_file) {
                    getline(_input_file, _line);
                    if (_input_file.eof()) break;

                    vector<string> _row;
                    boost::trim(_line);
                    boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);

                    int _i_index = boost::lexical_cast<int>(_row[0]);
                    int _j_index = boost::lexical_cast<int>(_row[1]);
                    //cout << "Vxc element [" << _i_index << ":" << _j_index << "] " << boost::lexical_cast<double>( _row[2] ) << endl;
                    _vxc(_i_index - 1, _j_index - 1) = boost::lexical_cast<double>(_row[2]);
                }

                LOG(logDEBUG, *_pLog) << "Done parsing" << flush;
                _input_file.close();
            }




            return true;
        }

        /**
         * Converts the Gaussian data stored in the Orbitals object to GW input format
         * This includes - writing out the orbitals with additional zeros and permutations
         *               - rewriting the Vxc matrix from cartesian to spherical AOs
         *               - forming the expectation value of the MOs with Vxc
         *               - writing out basis set information
         */
        bool Gaussian::ConvertToGW(Orbitals * _orbitals) {

            LOG(logDEBUG, *_pLog) << "Converting Gaussian to GW " << flush;
            // prepare file names
            string _orb_file_name_full = _run_dir + "/orbitals.votca";
            string _vxc_file_name_full = _run_dir + "/vxc.votca";

            // reload the basis set
            list<string> elements;
            BasisSet bs;
            // string _basisset_name("ubecppol");

            bs.LoadBasisSet(_basisset_name);
            LOG(logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;

            // rewriting the molecular orbitals
            ofstream _orb_file;
            _orb_file.open(_orb_file_name_full.c_str());
            _orb_file.precision(8);
            // getting the basis set sizes
            //int _basis_size      = _orbitals->getBasisSetSize();
            std::vector<double>::size_type _basis_size = _orbitals->getBasisSetSize();
            //int _cart_basis_size = _orbitals->_vxc.size1();
            //ub::matrix<double>::size_type _cart_basis_size = _orbitals->AOVxc().size1();
            //cout << "\nSpherical basis size is " << _basis_size << endl;
            //cout << "\nCartesian basis size is " << _cart_basis_size << endl;

            ub::matrix<double> mo_coefficients = (*_orbitals->getOrbitals());
            ub::vector<double> energies = (*_orbitals->getEnergies());

            // Sanity checks
            if (mo_coefficients.size1() != _basis_size) {
                cerr << "Incompatible basis size in molecular orbitals " << mo_coefficients.size1() << " vs basis set " << _basis_size;
                throw std::runtime_error("Conversion failed!");
                return false;
            }
            if (energies.size() != _basis_size) {
                cerr << "Incompatible basis size in energies " << energies.size() << " vs basis set " << _basis_size;
                throw std::runtime_error("Conversion failed!");
                return false;
            }

            // orbital energies in GW are in Rydberg
            energies = 2.0 * energies;

            // get atoms from orbitals and define iterator
            std::vector< QMAtom* > atoms = (*_orbitals->getAtoms());
            vector< QMAtom* >::iterator ita;
            LOG(logDEBUG, *_pLog) << "Rewriting molecular orbitals " << flush;
            // Loop over all molecular orbitals
            for (unsigned _i_orbital = 0; _i_orbital < _basis_size; _i_orbital++) {
                _orb_file << _i_orbital + 1 << " " << FortranFormat(energies(_i_orbital)) << endl;
                int _i_coef_qc = 0;
                int _i_coef_gw = 0;

                // Loop over all atoms
                for (ita = atoms.begin(); ita < atoms.end(); ita++) {

                    string element_name = (*ita)->type;
                    Element* element = bs.getElement(element_name);
                    //cout << "Atom " << element_name << endl;
                    // go through all shell types of this element
                    for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                        string shell_type = (*its)->getType();

                        // write out coefficients as needed
                        if (shell_type == "S") {
                            // write the s function coefficient
                            _orb_file << _i_orbital + 1 << " " << _i_coef_gw + 1 << " " << mo_coefficients(_i_orbital, _i_coef_qc) << endl;
                        } else if (shell_type == "P") {
                            // add one zero for unused s function
                            _orb_file << _i_orbital + 1 << " " << _i_coef_gw + 1 << " " << 0.0 << endl;
                            // write the px, py, pz function coefficients
                            for (int j = 0; j < 3; j++) {
                                _orb_file << _i_orbital + 1 << " " << _i_coef_gw + j + 2 << " " << mo_coefficients(_i_orbital, _i_coef_qc + j) << endl;
                            }
                        } else if (shell_type == "D") {
                            // add four zeros for unused s, px, py, pz functions
                            for (int j = 1; j < 5; j++) {
                                _orb_file << _i_orbital + 1 << " " << _i_coef_gw + j << " " << 0.0 << endl;
                            }
                            /* write the d function coefficients
                             * Gaussian has 5 spherical d functions (3z2-r2,xz,yz,x2-y2,xy)
                             * isogwa   has 5 spherical d functions (xz,yz,xy,3zz-rr,xx-yy)
                             */
                            _orb_file << _i_orbital + 1 << " " << _i_coef_gw + 5 << " " << mo_coefficients(_i_orbital, _i_coef_qc + 1) << endl;
                            _orb_file << _i_orbital + 1 << " " << _i_coef_gw + 6 << " " << mo_coefficients(_i_orbital, _i_coef_qc + 2) << endl;
                            _orb_file << _i_orbital + 1 << " " << _i_coef_gw + 7 << " " << mo_coefficients(_i_orbital, _i_coef_qc + 4) << endl;
                            _orb_file << _i_orbital + 1 << " " << _i_coef_gw + 8 << " " << mo_coefficients(_i_orbital, _i_coef_qc) << endl;
                            _orb_file << _i_orbital + 1 << " " << _i_coef_gw + 9 << " " << mo_coefficients(_i_orbital, _i_coef_qc + 3) << endl;

                            // add one zeros for unused s* function
                            _orb_file << _i_orbital + 1 << " " << _i_coef_gw + 10 << " " << 0.0 << endl;
                        } else {
                            cerr << "Conversion of shell type " << shell_type << " is not implemented!";
                            throw std::runtime_error("Conversion failed!");
                            return false;
                        }
                        _i_coef_gw += NumbfGW(shell_type);
                        _i_coef_qc += NumbfQC(shell_type);
                    }
                }
            }
            _orb_file.close();


            /* The AO Vxc matrix is expressed in cartesian Gaussians, while the MO are 
             * in spherical Gaussians. The AO matrix has first to be transformed to
             * sphericals.
             */

            ub::matrix<double> vxc = ub::zero_matrix<double> (_basis_size, _basis_size);
            const ub::matrix<double> vxc_cart = _orbitals->AOVxc();
            vector< QMAtom* >::iterator jta;

            int _isp = 0;
            int _jsp = 0;
            int _ica = 0;
            int _jca = 0;
            LOG(logDEBUG, *_pLog) << "Rewriting Vxc from cartesian to spherical functions " << flush;
            // loop over all i_atoms
            for (ita = atoms.begin(); ita < atoms.end(); ita++) {

                string i_element_name = (*ita)->type;
                Element* i_element = bs.getElement(i_element_name);
                // go through all shell types of this element
                for (Element::ShellIterator its = i_element->firstShell(); its != i_element->lastShell(); its++) {
                    string i_shell_type = (*its)->getType();
                    _jsp = 0;
                    _jca = 0;
                    // loop over all j_atoms
                    for (jta = atoms.begin(); jta < atoms.end(); jta++) {

                        string j_element_name = (*jta)->type;
                        Element* j_element = bs.getElement(j_element_name);
                        // go through all shell types of this element
                        for (Element::ShellIterator jts = j_element->firstShell(); jts != j_element->lastShell(); jts++) {
                            string j_shell_type = (*jts)->getType();

                            //cout << i_shell_type << "   "  << j_shell_type << endl;
                            if ((i_shell_type == "S" || i_shell_type == "P") && (j_shell_type == "S" || j_shell_type == "P")) {
                                // <s,p|vxc|s,p> elements are just copy and paste jobs
                                for (int i = 0; i < NumbfQC(i_shell_type); i++) {
                                    for (int j = 0; j < NumbfQC_cart(j_shell_type); j++) {
                                        vxc(_isp + i, _jsp + j) = vxc_cart(_ica + i, _jca + j);
                                    }
                                }
                            } else if ((i_shell_type == "S" || i_shell_type == "P") && j_shell_type == "D") {
                                // <s,p|V_xc|d> must be rewritten
                                for (int i = 0; i < NumbfQC(i_shell_type); i++) {
                                    // <s,p|Vxc|d3z2-r2>
                                    vxc(_isp + i, _jsp) = -0.5 * vxc_cart(_ica + i, _jca)
                                            - 0.5 * vxc_cart(_ica + i, _jca + 1)
                                            + vxc_cart(_ica + i, _jca + 2);
                                    // <s,p|Vxc|dxz>
                                    vxc(_isp + i, _jsp + 1) = vxc_cart(_ica + i, _jca + 4);

                                    // <s,p|Vxc|dyz>
                                    vxc(_isp + i, _jsp + 2) = vxc_cart(_ica + i, _jca + 5);

                                    // <s,p|Vxc|dx2-y2>
                                    vxc(_isp + i, _jsp + 3) = 0.5 * sqrt(3.0)
                                            *(vxc_cart(_ica + i, _jca)
                                            - vxc_cart(_ica + i, _jca + 1));
                                    // <s,p|Vxc|dxy>
                                    vxc(_isp + i, _jsp + 4) = vxc_cart(_ica + i, _jca + 3);
                                }
                            } else if (i_shell_type == "D" && (j_shell_type == "S" || j_shell_type == "P")) {
                                // <d|V_xc|s,p> must be rewritten
                                for (int j = 0; j < NumbfQC(j_shell_type); j++) {
                                    // <d3z2-r2|v_xc|s,p>
                                    vxc(_isp, _jsp + j) = -0.5 * vxc_cart(_ica, _jca + j)
                                            - 0.5 * vxc_cart(_ica + 1, _jca + j)
                                            + vxc_cart(_ica + 2, _jca + j);
                                    // <dxz|v_xc|s,p> 
                                    vxc(_isp + 1, _jsp + j) = vxc_cart(_ica + 4, _jca + j);

                                    // <dyz|v_xc|s,p>
                                    vxc(_isp + 2, _jsp + j) = vxc_cart(_ica + 5, _jca + j);

                                    // <dx2-y2|v_xc|s,p>
                                    vxc(_isp + 3, _jsp + j) = 0.5 * sqrt(3.0)
                                            *(vxc_cart(_ica, _jca + j)
                                            - vxc_cart(_ica + 1, _jca + j));
                                    // <dxy|v_xc|s,p>
                                    vxc(_isp + 4, _jsp + j) = vxc_cart(_ica + 3, _jca + j);

                                }

                            } else if (i_shell_type == "D" && j_shell_type == "D") {
                                // <d|Vxc|d> has to be rewritten

                                // <d3z2-r2|Vxc|d3z2-r2>
                                vxc(_isp, _jsp) = 0.25 * (vxc_cart(_ica, _jca)
                                        + vxc_cart(_ica, _jca + 1)
                                        + vxc_cart(_ica + 1, _jca)
                                        + vxc_cart(_ica + 1, _jca + 1))
                                        - 0.5 * (vxc_cart(_ica, _jca + 2)
                                        + vxc_cart(_ica + 1, _jca + 2)
                                        + vxc_cart(_ica + 2, _jca)
                                        + vxc_cart(_ica + 2, _jca + 1))
                                        + vxc_cart(_ica + 2, _jca + 2);

                                // <d3z2-r2|Vxc|dxz>
                                vxc(_isp, _jsp + 1) = -0.5 * (vxc_cart(_ica, _jca + 4)
                                        + vxc_cart(_ica + 1, _jca + 4))
                                        + vxc_cart(_ica + 2, _jca + 4);

                                // <d3z2-r2|Vxc|dyz>
                                vxc(_isp, _jsp + 2) = -0.5 * (vxc_cart(_ica, _jca + 5)
                                        + vxc_cart(_ica + 1, _jca + 5))
                                        + vxc_cart(_ica + 2, _jca + 5);

                                // <d3z2-r2|Vxc|dx2-y2>
                                vxc(_isp, _jsp + 3) = 0.25 * sqrt(3.0)
                                        *(-vxc_cart(_ica, _jca)
                                        + vxc_cart(_ica, _jca + 1)
                                        - vxc_cart(_ica + 1, _jca)
                                        + vxc_cart(_ica + 1, _jca + 1)
                                        + 2.0 * vxc_cart(_ica + 2, _jca)
                                        - 2.0 * vxc_cart(_ica + 2, _jca + 1));

                                //  <d3z2-r2|Vxc|dxy>  
                                vxc(_isp, _jsp + 4) = -0.5 * (vxc_cart(_ica, _jca + 3)
                                        + vxc_cart(_ica + 1, _jca + 3))
                                        + vxc_cart(_ica + 2, _jca + 3);

                                // <dxz|V|d3z2-r2>
                                vxc(_isp + 1, _jsp) = -0.5 * (vxc_cart(_ica + 4, _jca)
                                        + vxc_cart(_ica + 4, _jca + 1))
                                        + vxc_cart(_ica + 4, _jca + 2);

                                // <dxz|V|dxz>
                                vxc(_isp + 1, _jsp + 1) = vxc_cart(_ica + 4, _jca + 4);

                                // <dxz|V|dyz>
                                vxc(_isp + 1, _jsp + 2) = vxc_cart(_ica + 4, _jca + 5);

                                // <dxz|V|dx2-y2>
                                vxc(_isp + 1, _jsp + 3) = 0.5 * sqrt(3.0)
                                        *(vxc_cart(_ica + 4, _jca)
                                        - vxc_cart(_ica + 4, _jca + 1));

                                // <dxz|V|dxy>
                                vxc(_isp + 1, _jsp + 4) = vxc_cart(_ica + 4, _jca + 3);

                                // <dyz|V|d3z2-r2>
                                vxc(_isp + 2, _jsp) = -0.5 * (vxc_cart(_ica + 5, _jca)
                                        + vxc_cart(_ica + 5, _jca + 1))
                                        + vxc_cart(_ica + 5, _jca + 2);

                                // <dyz|V|dxz>
                                vxc(_isp + 2, _jsp + 1) = vxc_cart(_ica + 5, _jca + 4);

                                // <dyz|V|dyz>
                                vxc(_isp + 2, _jsp + 2) = vxc_cart(_ica + 5, _jca + 5);

                                // <dyz|V|dx2-y2>
                                vxc(_isp + 2, _jsp + 3) = 0.5 * sqrt(3.0)
                                        * (vxc_cart(_ica + 5, _jca)
                                        - vxc_cart(_ica + 5, _jca + 1));

                                // <dyz|V|dxy>
                                vxc(_isp + 2, _jsp + 4) = vxc_cart(_ica + 5, _jca + 3);

                                // <dx2-y2|V|d3z2-r2>
                                vxc(_isp + 3, _jsp) = 0.25 * sqrt(3.0)
                                        *(-vxc_cart(_ica, _jca)
                                        - vxc_cart(_ica, _jca + 1)
                                        + vxc_cart(_ica + 1, _jca)
                                        + vxc_cart(_ica + 1, _jca + 1)
                                        + 2.0 * vxc_cart(_ica, _jca + 2)
                                        - 2.0 * vxc_cart(_ica + 1, _jca + 2));

                                // <dx2-y2|V|dxz>
                                vxc(_isp + 3, _jsp + 1) = 0.5 * sqrt(3.0)
                                        *(vxc_cart(_ica, _jca + 4)
                                        - vxc_cart(_ica + 1, _jca + 4));
                                // <dx2-y2|V|dyz>
                                vxc(_isp + 3, _jsp + 2) = 0.5 * sqrt(3.0)
                                        *(vxc_cart(_ica, _jca + 5)
                                        - vxc_cart(_ica + 1, _jca + 5));

                                // <dx2-y2|V|dx2-y2>
                                vxc(_isp + 3, _jsp + 3) = 0.75 * (vxc_cart(_ica, _jca)
                                        - vxc_cart(_ica, _jca + 1)
                                        - vxc_cart(_ica + 1, _jca)
                                        + vxc_cart(_ica + 1, _jca + 1));

                                // <dx2-y2|V|dxy>
                                vxc(_isp + 3, _jsp + 4) = 0.5 * sqrt(3.0)
                                        *(vxc_cart(_ica, _jca + 3)
                                        - vxc_cart(_ica + 1, _jca + 3));

                                // <dxy|V|d3z2-r2>
                                vxc(_isp + 4, _jsp) = -0.5 * (vxc_cart(_ica + 3, _jca)
                                        + vxc_cart(_ica + 3, _jca + 1))
                                        + vxc_cart(_ica + 3, _jca + 2);
                                // <dxy|V|dxz>
                                vxc(_isp + 4, _jsp + 1) = vxc_cart(_ica + 3, _jca + 4);

                                // <dxy|V|dyz>
                                vxc(_isp + 4, _jsp + 2) = vxc_cart(_ica + 3, _jca + 5);


                                // <dxy|V|dx2-y2>
                                vxc(_isp + 4, _jsp + 3) = 0.5 * sqrt(3.0)
                                        *(vxc_cart(_ica + 3, _jca)
                                        - vxc_cart(_ica + 3, _jca + 1));
                                // <dxy|V|dxy>
                                vxc(_isp + 4, _jsp + 4) = vxc_cart(_ica + 3, _jca + 3);

                            }
                            // update matrix indices
                            _jsp += NumbfQC(j_shell_type);
                            _jca += NumbfQC_cart(j_shell_type);
                        } // j_shells
                    } // j_atoms
                    _isp += NumbfQC(i_shell_type);
                    _ica += NumbfQC_cart(i_shell_type);
                } // i_shells
            } // i_atoms

            /* Finally, with the rewritten AO Vxc matrix, we have to determine the
             * expectation values of the MOs with Vxc and write them to file
             */

            /* calculating the expectation values for all orbitals is overkill,
             * this should be limited to the range of orbitals included in the 
             * GW and BSE steps, according to options (LATER!)
             */
            LOG(logDEBUG, *_pLog) << "Calculating MO Vxc matrix " << flush;
            ub::matrix<double> vxc_expect = ub::zero_matrix<double> (_basis_size, _basis_size);
            ub::matrix<double> _temp = ub::zero_matrix<double> (_basis_size, _basis_size);
            _temp = ub::prod(vxc, ub::trans(mo_coefficients));
            vxc_expect = ub::prod(mo_coefficients, _temp);

            // write to file vxc.votca
            LOG(logDEBUG, *_pLog) << "Writing MO Vxc matrix " << flush;
            // output to file
            ofstream _vxc_file;
            _vxc_file.open(_vxc_file_name_full.c_str());
            for (unsigned _i_orbital = 0; _i_orbital < _basis_size; _i_orbital++) {
                for (unsigned _j_orbital = 0; _j_orbital < _basis_size; _j_orbital++) {
                    _vxc_file << _i_orbital + 1 << "  " << _j_orbital + 1 << "  " << FortranFormat(2.0 * vxc_expect(_i_orbital, _j_orbital)) << endl;
                }
            }
            _vxc_file.close();

            // 



            LOG(logDEBUG, *_pLog) << "Done converting to GW" << flush;
            return true;
        }

        string Gaussian::FortranFormat(const double &number) {
            stringstream _ssnumber;
            if (number >= 0) _ssnumber << " ";
            _ssnumber << setiosflags(ios::fixed) << setprecision(8) << std::scientific << number;
            std::string _snumber = _ssnumber.str();
            boost::replace_first(_snumber, "e", "D");
            return _snumber;
        }

        int Gaussian::NumbfGW(string shell_type) {
            int _nbf=0;
            if (shell_type == "S") {
                _nbf = 1;
            } else if (shell_type == "P") {
                _nbf = 4;
            } else if (shell_type == "D") {
                _nbf = 10;
            } else if (shell_type == "SP") {
                _nbf = 4;
            } else if (shell_type == "SPD") {
                _nbf = 10;
            }
            return _nbf;
        }

        int Gaussian::NumbfQC(string shell_type) {
            int _nbf=0;
            if (shell_type == "S") {
                _nbf = 1;
            } else if (shell_type == "P") {
                _nbf = 3;
            } else if (shell_type == "D") {
                _nbf = 5;
            } else if (shell_type == "SP") {
                _nbf = 4;
            } else if (shell_type == "SPD") {
                _nbf = 9;
            }
            return _nbf;
        }

        int Gaussian::NumbfQC_cart(string shell_type) {
            int _nbf=0;
            if (shell_type == "S") {
                _nbf = 1;
            } else if (shell_type == "P") {
                _nbf = 3;
            } else if (shell_type == "D") {
                _nbf = 6;
            } else if (shell_type == "SP") {
                _nbf = 4;
            } else if (shell_type == "SPD") {
                _nbf = 10;
            }
            return _nbf;
        }



    }
}
