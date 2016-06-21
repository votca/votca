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

            // NWChem file names
            string fileName = "system";

            _xyz_file_name = fileName + ".xyz";
            _input_file_name = fileName + ".inp";
            _log_file_name = fileName + ".log";
            _shell_file_name = fileName + ".sh";
            _orb_file_name = "WFNCOEF" ;

            string key = "package";
            string _name = options->get(key+".name").as<string> ();

            if ( _name != "cpmd" ) {
                cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error( "Wrong options file");
            }

            _executable =       options->get(key + ".executable").as<string> ();
            _charge =           options->get(key + ".charge").as<int> ();
            _spin =             options->get(key + ".spin").as<int> ();
            _options =          options->get(key + ".options").as<string> ();
            _memory =           options->get(key + ".memory").as<string> ();
            _threads =          options->get(key + ".threads").as<int> ();
            _scratch_dir =      options->get(key + ".scratch").as<string> ();
            _cleanup =          options->get(key + ".cleanup").as<string> ();





            //restart?
            if (options->exists(key + ".restart")) {
                _rsrt=true;
                _rsrt_kwds = options->get(key + ".restart").as<std::string> ();
            }
            else _rsrt=false;

            //optimize wavefunction?
            if (options->exists(key + ".optimizewf")) {
                _optWF=true;
                _convCutoff = options->get(key + ".optimizewf").as<double> ();
            }
            else _optWF=false;


            //functional
            if (options->exists(key + ".functional")) {
                _functional = options->get(key + ".functional").as<std::string> ();
            }
            else throw std::runtime_error("No functional specified");
            
            //symmetry
            _symmetry=0;
            if (options->exists(key + ".symmetry")) {
                _symmetry = options->get(key + ".symmetry").as<int> ();
            }
            else
                LOG(logDEBUG, *_pLog) << "CPMD: no symmetry provided, assuming simple cubic." << flush;
            
            //cell
            _cell="20.0   1.0   1.0  0.0  0.0  0.0";
            if (options->exists(key + ".cell")) {
                _cell = options->get(key + ".cell").as<std::string> ();
            }
            else
                LOG(logDEBUG, *_pLog) << "CPMD: no cell provided, assuming cube with side length of 20 Bohr." << flush;
            
            //plane wave cutoff
            _pwCutoff=80.0;
            if (options->exists(key + ".pwcutoff")) {
                _pwCutoff = options->get(key + ".pwcutoff").as<double> ();
            }
            else
                LOG(logDEBUG, *_pLog) << "CPMD: no plane wave cutoff provided, assuming "<< _pwCutoff <<" Ry." << flush;

            //output electrostatic potential?
            _elpot=false;
            if (options->exists(key + ".elpot")) {
                _elpot=options->get(key + ".elpot").as<bool> ();
            }

            //project wavefunction onto atomic orbitals?
            _projectWF=false;
            if (options->exists(key + ".projectwf")) {
                _projectWF=options->get(key + ".projectwf").as<bool> ();
            }

            //do population analysis? required to have access to WF coefs in atom-centric basis
            //requires _projectWF
            _popAnalysis=false;
            if (options->exists(key + ".popanalysis")) {
                _popAnalysis=options->get(key + ".popanalysis").as<bool> ();
            }

            //get density and overlap matrices during post processing?
            //requires _popAnalysis and _projectWF
            _getMat=false;
            if (options->exists(key + ".getmatrices")) {
                _getMat=options->get(key + ".getmatrices").as<bool> ();
            }

            if(_getMat) _popAnalysis=true;
            if(_popAnalysis) _projectWF=true;

            if(_projectWF && _optWF){
                cerr << "Error: Wavefunction optimization and projection onto atom-centric orbitals can not be done together.\nCPMD would crash.\n";
                cerr << "Do Wavefunction optimization first and then do projection/population analysis\n";
                cerr << "in a separate run with <restart>WAVEFUNCTION</restart>\n" << flush;
                LOG(logDEBUG, *_pLog) << "CPMD: Wavefunction optimization and projection onto atom-centric orbitals can not be done together." << flush;
                throw std::runtime_error("Mutually exclusive options");
            }

        }

        bool Cpmd::WriteInputFile(std::vector<Segment* > segments, Orbitals* orbitals_guess) {

            std::vector< Atom* > _atoms;
            std::vector< Atom* > ::iterator ait;
            std::vector< Segment* >::iterator sit;
            std::string temp_suffix = "/id";
            std::string scratch_dir_backup = _scratch_dir;

            ofstream _com_file;

            std::string _com_file_name_full = _run_dir + "/" + _input_file_name;

            _com_file.open(_com_file_name_full.c_str());

            // header
            _com_file << "&INFO\nGenerated by VOTCA\n&END" << endl;


            //control
            _com_file << "\n&CPMD\n";
            if(_rsrt) _com_file << "  RESTART " << _rsrt_kwds << endl;  //restart
            if(_optWF){                                                 //optimize WF
                _com_file << "  OPTIMIZE WAVEFUNCTION" << endl;
                _com_file << "  CONVERGENCE ORBITALS" << endl;
                _com_file << "  " << FortranFormat(_convCutoff) << endl;
            }
            if(_elpot){                                                 //output electrostatic potential
                _com_file << "  ELECTROSTATIC POTENTIAL" << endl;
                _com_file << "  RHOOUT" << endl;
            }
            if(_projectWF)
                _com_file << "  PROPERTIES" << endl;
            _com_file << "&END" << endl;

            //functional
            _com_file << "\n&DFT\n";
            _com_file << "  FUNCTIONAL " << _functional << endl;
            _com_file << "  GC-CUTOFF" << endl;
            _com_file << "   1.0d-06" << endl;
            _com_file << "&END" << endl;
            
            //cell
            _com_file << "\n&SYSTEM\n";
            _com_file << "  SYMMETRY" << endl;
            _com_file << "   " << _symmetry <<endl;
            _com_file << "  CELL" << endl;
            _com_file << "   " << _cell <<endl;
            _com_file << "  CUTOFF" << endl;
            _com_file << "   " << FortranFormat(_pwCutoff) <<endl;
            _com_file << "&END" << endl;
            
            //properties
            if(_projectWF){
                _com_file << "\n&PROP\n";
                _com_file << "  PROJECT WAVEFUNCTION" << endl;
                if(_popAnalysis){
                    _com_file << "  CHARGES" << endl;
                    _com_file << "  POPULATION ANALYSIS MULLIKEN" << endl;
                }
                _com_file << "&END" << endl;
            }
            
            //basis
            _com_file << "\n&BASIS\n";
            WriteBasisSet(segments, _com_file);
            _com_file << "&END" << endl;
            
            //atoms
            _com_file << "\n&ATOMS\n";
            #warning "TODO: put atomic coordinates into the input file"
            _com_file << "&END" << endl;
            
            
            
            _com_file << endl;
            _com_file.close();

            return true;
        }
        
        
        /**
         * Writes the basis set files to disk in a format that CPMD can understand
         */
        void Cpmd::WriteBasisSet(std::vector<Segment* > segments, ofstream &_com_file) {
            
            std::vector< Atom* > _atoms;
            std::vector< Atom* > ::iterator ait;
            std::vector< Segment* >::iterator sit;
            list<std::string> elements;
            BasisSet bs;

            bs.LoadBasisSet(_basisset_name);
            LOG(logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;

            for (sit = segments.begin(); sit != segments.end(); ++sit) {

                std::vector< Atom* > atoms = (*sit)-> Atoms();
                std::vector< Atom* >::iterator it;

                for (it = atoms.begin(); it < atoms.end(); it++) {

                    std::string element_name = (*it)->getElement();

                    list<std::string>::iterator ite;
                    ite = find(elements.begin(), elements.end(), element_name);

                    if (ite == elements.end()) {
                        elements.push_back(element_name);

                        Element* element = bs.getElement(element_name);
                        
                        std::string _short_el_file_name = element_name + "_" + _basisset_name + ".basis";
                        std::string _el_file_name = _run_dir + "/" + _short_el_file_name;
                        
                        
                        //write the element to the input file
                        _com_file << "*" << _short_el_file_name << " " << std::distance(element->firstShell(), element->lastShell()) << " GAUSSIAN"<<endl;
                        _com_file << "   ";
                        
                        //create the .basis file
                        ofstream _el_file;
                        _el_file.open(_el_file_name.c_str());
                        
                        //comment
                        _el_file << element_name << " with the "<< _basisset_name << " basis." << endl;
                        
                        //Lmax
                        _el_file << element->getLmax() << endl;
                        
                        //sort shells by L
                        for (int L=0; L <= element->getLmax(); L++)
                        {
                            _el_file << "  Functions for l="<<L<<endl;
                            
                            std::vector<Shell*> Lshells;
                            
                            int ndecays=0;
                            for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                Shell* shell = (*its);
                                int Ls=shell->getLmax();
                                if(shell->getType().size()>1){
                                    cerr << "CPMD does not support " << shell->getType() << "basis functions." << endl;
                                    cerr << "Please break the basis set into basis functions with only one L-value each." << endl << flush;
                                    LOG(logDEBUG, *_pLog) << "CPMD: multi-L basis functions not supported." << flush;
                                    throw std::runtime_error("Unsupported basis function");
                                }

                                //For now assume all shells have only one L-value.
                                //Can decompose the basis set into such shells later, in another place.
                                //Any subsequent analysis will have to use the decomposed basis set too.
                                if (Ls==L) //this shell has the correct L
                                {
                                    ndecays+=shell->getSize();
                                    Lshells.push_back(shell);
                                    
                                    //write the shell's L-value to the input file
                                    _com_file << L << " ";
                                }
                            }
                            
                            _el_file << "  " << Lshells.size()<< " " << ndecays << endl;
                            _el_file << endl;
                            
                            //decays
                            ios::fmtflags old_settings = _el_file.flags();
                            _el_file << std::scientific << std::setprecision(6);
                            _el_file << "  ";
                            for (Element::ShellIterator its = Lshells.begin(); its !=  Lshells.end(); its++)
                            {
                                Shell* shell = (*its);
                                for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                    GaussianPrimitive* gaussian = *itg;
                                    _el_file << gaussian->decay << "\t";
                                }
                            }
                            _el_file << endl;
                            
                            //coefficients (scale*contraction)
                            int gs=0; //number of decays already handled
                            for (Element::ShellIterator its = Lshells.begin(); its !=  Lshells.end(); its++)
                            {
                                Shell* shell = (*its);
                                if(shell->getSize()!=0) //there are gaussians in this shell
                                {
                                    int gi=0; //index of the current decay
                                    _el_file << "  ";
                                    //output zeros for all decays already handled
                                    for(gi=0; gi<gs; gi++)
                                    {
                                        _el_file << 0.0 << "\t";
                                    }
                                    //output coefficients for this shell's gaussians
                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                        _el_file << shell->getScale() * gaussian->contraction[0] << "\t";
                                        gi++;
                                    }
                                    gs+=shell->getSize();
                                    //output zeros till the end of decays
                                    for(;gi<ndecays; gi++)
                                    {
                                        _el_file << 0.0 << "\t";
                                    }
                                    _el_file << endl;
                                }
                            }
                            _el_file.flags(old_settings);
                        }

                        _el_file << endl;
                        _el_file.close();
                        
                        _com_file<<endl;
                    }
                }
            }
        }
        

        /**
         * Runs the CPMD job.
         */
        bool Cpmd::Run() {

            LOG(logDEBUG, *_pLog) << "CPMD: running [" << _executable << " " << _input_file_name << "]" << flush;

            if (std::system(NULL)) {
                std::string _command;
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

        
        std::string Cpmd::FortranFormat(const double &number) {
            std::stringstream _ssnumber;
            if (number >= 0) _ssnumber << " ";
            _ssnumber << setiosflags(ios::fixed) << setprecision(8) << std::scientific << number;
            std::string _snumber = _ssnumber.str();
            boost::replace_first(_snumber, "e", "D");
            return _snumber;
        }



    }
}
