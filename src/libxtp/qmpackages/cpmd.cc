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

#include "cpmd.h"
#include <votca/ctp/segment.h>

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


            //functional and pseudopotentials
            if (options->exists(key + ".functional")) {
                if (!(options->exists(key + ".functional"))) throw std::runtime_error("Functional name missing");
                _functional = options->get(key + ".functional.name").as<std::string> ();
                
                list<Property*>::iterator pit;
                std::list<Property *> props = options->Select(key + ".functional.pseudopotentials.*");
                for (pit = props.begin(); pit != props.end(); ++pit)
                {
                    Property* p= *pit;
                    _ppFileNames[p->name()]=p->as<std::string> ();
                }
                props = options->Select(key + ".functional.l.*");
                for (pit = props.begin(); pit != props.end(); ++pit)
                {
                    Property* p= *pit;
                    _ppLData[p->name()]=p->as<std::string> ();
                }
                
            }
            else throw std::runtime_error("No functional and pseudopotentials specified");
            
            //symmetry
            _symmetry=0;
            if (options->exists(key + ".symmetry")) {
                _symmetry = options->get(key + ".symmetry").as<int> ();
            }
            else
                LOG(ctp::logDEBUG, *_pLog) << "CPMD: no symmetry provided, assuming simple cubic." << flush;
            
            //cell
            _cell="20.0   1.0   1.0  0.0  0.0  0.0";
            if (options->exists(key + ".cell")) {
                _cell = options->get(key + ".cell").as<std::string> ();
            }
            else
                LOG(ctp::logDEBUG, *_pLog) << "CPMD: no cell provided, assuming cube with side length of 20 Bohr." << flush;
            
            //plane wave cutoff
            _pwCutoff=80.0;
            if (options->exists(key + ".pwcutoff")) {
                _pwCutoff = options->get(key + ".pwcutoff").as<double> ();
            }
            else
                LOG(ctp::logDEBUG, *_pLog) << "CPMD: no plane wave cutoff provided, assuming "<< _pwCutoff <<" Ry." << flush;

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
            
            //basis set file name
            if (options->exists(key + ".basisset")) {
            _basisset_name=options->get(key + ".basisset").as<string> ();
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
                LOG(ctp::logDEBUG, *_pLog) << "CPMD: Wavefunction optimization and projection onto atom-centric orbitals can not be done together." << flush;
                throw std::runtime_error("Mutually exclusive options");
            }

        }

        bool Cpmd::WriteInputFile(std::vector<ctp::Segment* > segments, Orbitals* orbitals_guess) {

            std::vector< ctp::Atom* > _atoms;
            std::vector< ctp::Atom* > ::iterator ait;
            std::vector< ctp::Segment* >::iterator sit;
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
            if(_projectWF){
                _com_file << "\n&BASIS\n";
                WriteBasisSet(segments, _com_file);
                _com_file << "&END" << endl;
            }
            
            //atoms
            _com_file << "\n&ATOMS\n";
                //find how many atoms of each element there are
            //list<std::string> elements;
            for (sit = segments.begin(); sit != segments.end(); ++sit) {

                _atoms = (*sit)-> Atoms();

                for (ait = _atoms.begin(); ait < _atoms.end(); ait++) {

                    std::string element_name = (*ait)->getElement();
                    list<std::string>::iterator ite;
                    ite = find(_elements.begin(), _elements.end(), element_name);
                    if (ite == _elements.end()) {            //this is the first atom of this element encountered
                        _elements.push_back(element_name);
                        _nAtomsOfElement[element_name]=1;
                    }
                    else
                    {
                        _nAtomsOfElement[element_name]++;
                    }
                }
            }
                //now loop over elements and store all atoms of that element
            list<std::string>::iterator ite;
            for (ite = _elements.begin(); ite != _elements.end(); ite++) {
                if(_ppFileNames.find(*ite)==_ppFileNames.end()) {
                    cerr << "Error: Element "<<(*ite)<<" has no pseudopotential specified in CPMD options file.\n" << flush;
                    throw std::runtime_error("Encountered element with no pseudopotential.\n");
                }
                else{
                    _com_file << "*" << _ppFileNames[(*ite)] << endl; //store name of the pseudopotential file and it's read options
                    if(_ppLData.find(*ite)==_ppLData.end()) {
                        cerr << "Warning: Element "<<(*ite)<<" has no angular momentum data (<l></l>) specified in CPMD options file.\n\tAttempting to read it from basis set. This may produce errors.\n" << flush;
                        if(_basisset_name.empty()){
                            cerr << "Error: Basis set file not specified.\n" << flush;
                            throw std::runtime_error("Encountered element with no angular momentum data.\n");
                        }
                        else{
                            BasisSet _bs;
                            _bs.LoadBasisSet(_basisset_name);
                            int Lmax = 0;
                            //find Lmax by checking all shells of the element
                            Element* el=_bs.getElement(*ite);
                            for (Element::ShellIterator its = el->firstShell(); its != el->lastShell(); its++) {
                                        int Ls=(*its)->getLmax();
                                        if(Lmax<Ls) Lmax=Ls;
                            }
                            _com_file << "   "<<Lmax<<" "<<Lmax<<" "<<Lmax<< endl; //LMAX LOC SKIP
                        }
                    }
                    else{
                        _com_file << "   "<<_ppLData[(*ite)]<< endl; //LMAX LOC SKIP
                    }
                    _com_file << "   "<< _nAtomsOfElement[(*ite)] <<endl;  //# atoms of element
                    
                    //store atomic positions
                    for (sit = segments.begin(); sit != segments.end(); ++sit) {
                        _atoms = (*sit)-> Atoms();
                        for (ait = _atoms.begin(); ait < _atoms.end(); ait++) {
                            if((*ait)->getElement().compare(*ite)==0){     //this element
                                vec pos = (*ait)->getQMPos();
                                _com_file << "   ";
                                _com_file << setw(12) << setiosflags(ios::fixed) << setprecision(5) << conv::nm2bohr*pos.getX() << "   ";
                                _com_file << setw(12) << setiosflags(ios::fixed) << setprecision(5) << conv::nm2bohr*pos.getY() << "   ";
                                _com_file << setw(12) << setiosflags(ios::fixed) << setprecision(5) << conv::nm2bohr*pos.getZ() << "   ";
                                _com_file << endl;
                            }
                        }
                    }
                }
            }
            
            //#warning "TODO: copy pseudopotentials to the _run_dir"
            _com_file << "&END" << endl;
            
            
            
            _com_file << endl;
            _com_file.close();

            return true;
        }
        
        
        /**
         * Writes the basis set files to disk in a format that CPMD can understand
         */
        void Cpmd::WriteBasisSet(std::vector<ctp::Segment* > segments, ofstream &_com_file) {
            
            std::vector< ctp::Atom* > _atoms;
            std::vector< ctp::Atom* > ::iterator ait;
            std::vector< ctp::Segment* >::iterator sit;
            list<std::string> elements;
            
            BasisSet _bs;
            _bs.LoadBasisSet(_basisset_name);
            LOG(ctp::logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;

            for (sit = segments.begin(); sit != segments.end(); ++sit) {

                std::vector< ctp::Atom* > atoms = (*sit)-> Atoms();
                std::vector< ctp::Atom* >::iterator it;

                for (it = atoms.begin(); it < atoms.end(); it++) {

                    std::string element_name = (*it)->getElement();

                    list<std::string>::iterator ite;
                    ite = find(elements.begin(), elements.end(), element_name);

                    if (ite == elements.end()) {
                        elements.push_back(element_name);

                        Element* element = _bs.getElement(element_name);
                        
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
                        int Lmax = 0;
                        //find Lmax by checking all shells of the element
                        for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                    int Ls=(*its)->getLmax();
                                    if(Lmax<Ls) Lmax=Ls;
                        }
                        _el_file << Lmax+1 << endl; //number of L-values in this .basis file
                        
                        
                        //sort shells by L
                        for (int L=0; L <= element->getLmax(); L++)
                        {
                            std::vector<Shell*> Lshells;
                            
                            int ndecays=0;
                            for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                Shell* shell = (*its);
                                int Ls=shell->getLmax();
                                if(shell->getType().size()>1){
                                    cerr << "CPMD does not support " << shell->getType() << " basis functions." << endl;
                                    cerr << "Please break the basis set into basis functions with only one L-value each." << endl << flush;
                                    LOG(ctp::logDEBUG, *_pLog) << "CPMD: multi-L basis functions not supported." << flush;
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
                            
                            if(!Lshells.empty()){   //only write things if there are shells with this L
                                _el_file << "  Functions for l="<<L<<endl;
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
                                            _el_file << shell->getScale() * gaussian->contraction[L] << "\t";
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

            LOG(ctp::logDEBUG, *_pLog) << "CPMD: running [" << _executable << " " << _input_file_name << "]" << flush;

            if (std::system(NULL)) {
                std::string _command;
                _command = "cd " + _run_dir + "; rm -f LocalError*.log; " + _executable + " " + _input_file_name + ">" + _log_file_name;
                int check=std::system(_command.c_str());
                if (check==-1){
                    LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                    return false;
                }
                if (CheckLogFile()) {
                    LOG(ctp::logDEBUG, *_pLog) << "CPMD: finished job" << flush;
                    return true;
                } else {
                    LOG(ctp::logDEBUG, *_pLog) << "CPMD: job failed" << flush;
                }
            } else {
                LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
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
           

            std::string _full_name = (arg_path / _run_dir / _log_file_name).c_str();
            ifstream _input_file(_full_name.c_str());

            if (_input_file.fail()) {
                LOG(ctp::logERROR, *_pLog) << "CPMD: " << _full_name << " is not found." << endl << flush;
                return false;
            };

            //Use brute force. Search every line for the termination string.
            //It doesn't appear at the very end, like in gaussian
            std::string::size_type self_energy_pos=std::string::npos;
            std::string _line;
            do {
                getline(_input_file, _line);
                self_energy_pos=_line.find("PROGRAM CPMD ENDED AT");
            } while (self_energy_pos==std::string::npos && !(_input_file.eof()));

            _input_file.close();

            if (self_energy_pos == std::string::npos) {
                LOG(ctp::logERROR, *_pLog) << "CPMD: " << _full_name << " is incomplete."<< endl << flush;
                return false;
            } else {
                LOG(ctp::logDEBUG,*_pLog) << "CPMD LOG is complete." <<endl << flush;
                return true;
            }
        }

        /**
         * Parses the CPMD Log file and stores data in the Orbitals object
         */
        bool Cpmd::ParseLogFile(Orbitals * _orbitals) {
            std::string _line;
            std::vector<std::string> results;
            
            LOG(ctp::logDEBUG, *_pLog) << "CPMD: parsing " << _log_file_name << flush;
            
            std::string _log_file_name_full = _log_file_name;
            if (_run_dir != "") _log_file_name_full = _run_dir + "/" + _log_file_name;
            
            // check if LOG file is complete
            if (!CheckLogFile()) return false;

            // save qmpackage name
            _orbitals->setQMpackage("cpmd");
            
            ifstream _input_file(_log_file_name_full.c_str());
            while (_input_file) {
                
                getline(_input_file, _line);
                boost::trim(_line);
                
                
                /*
                 * number of electrons
                 */
                std::string::size_type electrons_pos = _line.find("alpha electrons");
                if (electrons_pos != std::string::npos) {
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int _number_of_electrons = (int) boost::lexical_cast<double>(results.back());
                    _orbitals->setNumberOfElectrons(_number_of_electrons);
                    LOG(ctp::logDEBUG, *_pLog) << "Alpha electrons: " << _number_of_electrons << flush;
                }
                
            }
            LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
            _input_file.close();
            
            
            if(_projectWF){
                //MO coefficient and overlap matrices
                if(!loadMatrices(_orbitals)) return false;

                //atom info
                if(!(_orbitals->hasQMAtoms())){ //no atoms defined for the orbitals
                    //lets fill them in, in CPMD's order

                    //iterate over elements
                    list<std::string>::iterator ite;
                    int i=0;
                    for (ite = _elements.begin(); ite != _elements.end(); ite++, i++) {
                        for(int a=0; a<_NA[i]; a++){
                            _orbitals->AddAtom(*ite, 0, 0, 0, _ZV[i]); //store core charge in the atomic charge field
                        }
                    }
                }
                else
                {
                    cerr << "CPMD: _orbitals already has some atoms. Need to implement atom reordering for this case." << flush;
                    LOG(ctp::logDEBUG, *_pLog) << "CPMD: _orbitals already has some atoms. Need to implement atom reordering for this case." << flush;
                    throw std::runtime_error("Unimplemented case");
                    return false;
                }
            }
            
            
            //basis set
            if(_orbitals->hasGWbasis()){
                if(_orbitals->getGWbasis().compare(_basisset_name)!=0){
                    cerr << "CPMD: _orbitals already has a basis set and it does not match the basis set CPMD was initialized with." << flush;
                    LOG(ctp::logDEBUG, *_pLog) << "CPMD: _orbitals already has a basis set and it does not match the basis set CPMD was initialized with." << flush;
                    throw std::runtime_error("Basis set mismatch");
                    return false;
                }
            }
            else{
                _orbitals->setGWbasis(_basisset_name);
            }
            
            return true;

        }
        
        bool Cpmd::loadMatrices(Orbitals * _orbitals)
        {
            //check if WFNCOEF exists
            boost::filesystem::path arg_path;
            std::string _full_name = (arg_path / _run_dir / "WFNCOEF").c_str();
            ifstream wf_file(_full_name.c_str());
            if(wf_file.fail())
            {
                LOG(ctp::logERROR, *_pLog) << "CPMD: " << _full_name << " is not found." << endl << flush;
            }
            LOG(ctp::logDEBUG, *_pLog) << "CPMD: parsing " << _full_name << flush;
            
            //read WFNCOEF
            //variable names in all CAPS are variables from CPMD source code
            int count=0, endcount=0;
            wf_file.read((char*)&count, 4); //bytes in this record
            int bl=count;

            int NATTOT=0;
            wf_file.read((char*)&NATTOT, 4);    //number of basis functions
            bl-=4;
            _orbitals->setBasisSetSize(NATTOT);
            LOG(ctp::logDEBUG, *_pLog) << "Basis functions: " << NATTOT << flush;

            _NSP=0;                  //number of atom types
            wf_file.read((char*)&_NSP, 4);
            bl-=4;

            //double ZV[NSP];             //core charge
            //int NA[NSP], NUMAOR[NSP];
            if(_ZV!=NULL)
            {
                delete[] _ZV;
                delete[] _NA;
                delete[] _NUMAOR;
            }
            _ZV=new double[_NSP];
            _NA=new int[_NSP];
            _NUMAOR=new int[_NSP];
            if(_ZV==NULL || _NA==NULL || _NUMAOR==NULL)
                throw std::runtime_error("Memory allocation failed");
            for(int i=0; i<_NSP; i++)
            {
                wf_file.read((char*)&_ZV[i], 8);     //core charge of atom type i
                wf_file.read((char*)&_NA[i], 4);     //number of atoms of type i
                wf_file.read((char*)&_NUMAOR[i], 4); //number of atomic orbitals of atom type i
                bl-=8+4*2;
            }
            
            //check footer
            wf_file.read((char*)&endcount, 4);
            if(bl!=0 || endcount!=count){ //number of bytes read was wrong
                cerr << "CPMD: " << "could not parse record in "<< _full_name << endl << flush;
                LOG(ctp::logERROR, *_pLog) << "CPMD: " << "could not parse record in "<< _full_name << flush;
                throw std::runtime_error("IO error");
                return false;
            }
                
            //NEW RECORD
            wf_file.read((char*)&count, 4); //bytes in this record
            bl=count;

            int NUMORB=count/8/NATTOT;          //number of MOs (energy levels))
            LOG(ctp::logDEBUG, *_pLog) << "number of energy levels: " << NATTOT << flush;
                //resize the coefficient matrix
            ub::matrix<double> &mo_coefficients = _orbitals->MOCoefficients();
            mo_coefficients.resize(NUMORB, NATTOT);
            //double XXMAT[NATTOT][NUMORB];
            double XXMAT;
            for(int i=0; i<NUMORB; i++){
                for(int j=0; j<NATTOT; j++){  
                    wf_file.read((char*)&XXMAT, 8);
                    mo_coefficients(i,j)=XXMAT;
                    bl-=8;
                }
            }
            
            
            //check footer
            wf_file.read((char*)&endcount, 4);
            if(bl!=0 || endcount!=count){ //number of bytes read was wrong
                cerr << "CPMD: " << "could not parse record in "<< _full_name << endl << flush;
                LOG(ctp::logERROR, *_pLog) << "CPMD: " << "could not parse record in "<< _full_name << endl << flush;
                throw std::runtime_error("IO error");
                return false;
            }
            
            wf_file.close();
            LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
            
            
            
            
            //check if OVERLAP exists
            _full_name = (arg_path / _run_dir / "OVERLAP").c_str();
            ifstream ov_file(_full_name.c_str());
            if(ov_file.fail())
            {
                LOG(ctp::logERROR, *_pLog) << "CPMD: " << _full_name << " is not found." << endl << flush;
            }
            LOG(ctp::logDEBUG, *_pLog) << "CPMD: parsing " << _full_name << flush;
            
            //read OVERLAP
            count=0, endcount=0;
            ov_file.read((char*)&count, 4); //bytes in this record
            bl=count;
            
            if(NATTOT*NATTOT!=count/8)
            {
                cerr << "CPMD: " << "Number of basis functions in the overlap and coefficient matrices do not match."<< endl << flush;
                LOG(ctp::logERROR, *_pLog) << "CPMD: " << "Number of basis functions in the overlap and coefficient matrices do not match."<< endl << flush;
                throw std::runtime_error("IO error");
                return false;
            }
            
            
                //resize the overlap matrix
            ub::symmetric_matrix<double> &overlap = _orbitals->AOOverlap();
            overlap.resize(NATTOT);
            
            //read
            double XSMAT;
            for(int i=0; i<NATTOT; i++){
                for(int j=0; j<NATTOT; j++){  
                    ov_file.read((char*)&XSMAT, 8);
                    overlap(i,j)=XSMAT;
                    bl-=8;
                }
            }
            
            //check footer
            ov_file.read((char*)&endcount, 4);
            if(bl!=0 || endcount!=count){ //number of bytes read was wrong
                cerr << "CPMD: " << "could not parse record in "<< _full_name << endl << flush;
                LOG(ctp::logERROR, *_pLog) << "CPMD: " << "could not parse record in "<< _full_name << flush;
                throw std::runtime_error("IO error");
                return false;
            }
            
            ov_file.close();
            LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
            
            
           // #warning "TODO: Reorder Overlap and MOcoeefficient matrices so that the order of basis functions (like d_xx) matches the VOTCA order."
            
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
