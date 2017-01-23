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

#include "orca.h"
#include <votca/ctp/segment.h>
#include <votca/xtp/elements.h>
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



namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    namespace CTP = votca::ctp;
/******** Defining a class for l->number and number->l*****************************/
    class Angularmom
{
public:   
    Angularmom() { FillMaps(); };
   ~Angularmom() { };

    const int        &getLNum(std::string lname) const {return _LNum.at(lname); }
    const std::string     &getLName(int lnum) const {return _LName.at(lnum); }
    private:
    std::map<std::string, int> _LNum;
    std::map<int, std::string> _LName;
        inline void FillMaps(){
        
        FillLNum();
        FillLName();
        }
        
        inline void FillLName(){
    
        _LName[0]  = "S";
        _LName[1]  = "P";
        _LName[2]  = "D";
        _LName[3]  = "F";
        };
        
        inline void FillLNum(){
    
        _LNum["S"]  = 0;
        _LNum["P"]  = 1;
        _LNum["D"]  = 2;
        _LNum["F"]  = 3;
        };
     
 };
    
    
    /************************************************************/
void Orca::Initialize( Property *options ) {
    
    //good luck

     // Orca file names
    std::string fileName = "system";

    _xyz_file_name = fileName + ".xyz";
    _input_file_name = fileName + ".inp";
    _log_file_name = fileName + ".log"; 
    _shell_file_name = fileName + ".sh"; 
    _orb_file_name=fileName+".gbw";

    std::string key = "package";
    std::string _name = options->get(key+".name").as<std::string> ();
    
    if ( _name != "orca" ) {
        cerr << "Tried to use " << _name << " package. ";
        throw std::runtime_error( "Wrong options file");
    }
    
    _executable =       options->get(key + ".executable").as<std::string> ();
    _charge =           options->get(key + ".charge").as<int> ();
    _spin =             options->get(key + ".spin").as<int> ();
    _options =          options->get(key + ".options").as<std::string> ();
    _memory =           options->get(key + ".memory").as<std::string> ();
    _threads =          options->get(key + ".threads").as<int> ();
    _scratch_dir =      options->get(key + ".scratch").as<std::string> ();
    _basisset_name =    options->get(key + ".basisset").as<std::string> ();
    _cleanup =          options->get(key + ".cleanup").as<std::string> ();

    
    
    if (options->exists(key + ".outputVxc")) {
                _output_Vxc = options->get(key + ".outputVxc").as<bool> ();   
            }
             else _output_Vxc=false;
    if (_output_Vxc){
        throw std::runtime_error( "Sorry "+_name+" does not support Vxc output");
    }
  

    _write_basis_set=   options->get(key + ".writebasisset").as<bool> ();
    _write_pseudopotentials= options->get(key + ".writepseudopotentials").as<bool> ();
    

    // check if the optimize keyword is present, if yes, read updated coords
    std::string::size_type iop_pos = _options.find(" Opt"); /*optimization word in orca*/
    if (iop_pos != std::string::npos) {
        _is_optimization = true;
    } else
    {
        _is_optimization = false;
    }

    // check if the esp keyword is present, if yes, get the charges and save them
    iop_pos = _options.find(" chelpg");  /*electrostatic potential related to partial atomic charges I guess is chelpg in orca but check */
    if (iop_pos != std::string::npos) {
        _get_charges = true;
    } else
    {
        _get_charges = false;
    }

    // check if the charge keyword is present, if yes, get the self energy and save it
    iop_pos = _options.find("pointcharges");  /*??*/
    if (iop_pos != std::string::npos) {
        _get_self_energy = true;
        _write_charges = true;
    } else
    {
        _get_self_energy = false;
        _write_charges = false;
    }
    
    // check if the guess should be prepared, if yes, append the guess later
    _write_guess = false;
    iop_pos = _options.find("iterations 1 ");
    if (iop_pos != std::string::npos) _write_guess = true;
    iop_pos = _options.find("iterations 1\n");
    if (iop_pos != std::string::npos) _write_guess = true;
}    

/**
 * Prepares the *.inp file from a vector of segments
 * Appends a guess constructed from monomer orbitals if supplied
 */
bool Orca::WriteInputFile( std::vector<CTP::Segment* > segments, Orbitals* orbitals_guess )
{
    std::vector< CTP::Atom* > _atoms;
    std::vector< CTP::Atom* > ::iterator ait;
    std::vector< CTP::Segment* >::iterator sit;
    std::vector<std::string> results;
   // int qmatoms = 0;
    std::string temp_suffix = "/id";
    std::string scratch_dir_backup = _scratch_dir;
    std::ofstream _com_file;
    std::ofstream _crg_file;

    std::string _com_file_name_full = _run_dir + "/" + _input_file_name;
    std::string _crg_file_name_full = _run_dir + "/background.crg" ;
    

    _com_file.open ( _com_file_name_full.c_str() );
    // header 
   //_com_file << "* xyzfile  "  <<  _charge << " " << _spin << " " << _xyz_file_name << "\n" << endl;
    
    _com_file << "* xyz  "  <<  _charge << " " << _spin << endl;
    if (!_write_charges) {
                for (sit = segments.begin(); sit != segments.end(); ++sit) {
                    temp_suffix = temp_suffix + "_" + boost::lexical_cast<std::string>((*sit)->getId());
                    _atoms = (*sit)-> Atoms();

                    for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                        if ((*ait)->HasQMPart() == false) {
                            continue;
                        }

                        vec pos = (*ait)->getQMPos();
                        std::string name = (*ait)->getElement();

                        //fprintf(out, "%2s %4.7f %4.7f %4.7f \n"
                        _com_file << setw(3) << name.c_str()
                                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX()*10
                                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY()*10
                                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ()*10
                                << endl;
                    }
                }
                    _com_file << "* \n" << endl;
                _com_file << "%pal\n "  <<  "nprocs " <<  _threads  << "\nend" << "\n" << endl;
                
                // write basis set info to file
                if ( _write_basis_set) {
                     Elements _elements;
   
                    list<std::string> elements;
                    BasisSet bs;
                    bs.LoadBasisSet(_basisset_name);
                    LOG(CTP::logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;
                    ofstream _el_file;
                    std::string _el_file_name = _run_dir + "/" +  "system.bas";
                    _el_file.open(_el_file_name.c_str());
                    //_com_file << "@" << "system.bas" << endl;
                    _el_file << "$DATA" << endl;
                                                    
                    for (sit = segments.begin(); sit != segments.end(); ++sit) {
                        std::vector< CTP::Atom* > atoms = (*sit)-> Atoms();
                        std::vector< CTP::Atom* >::iterator it;
                        for (it = atoms.begin(); it < atoms.end(); it++) {
                            std::string element_name = (*it)->getElement();
                            list<std::string>::iterator ite;
                            ite = find(elements.begin(), elements.end(), element_name);
                            if (ite == elements.end()) {
                                elements.push_back(element_name);
                                Element* element = bs.getElement(element_name);
                                _el_file << _elements.getEleFull(element_name) << endl; 
                                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                    Shell* shell = (*its);
                                    _el_file  << shell->getType() << " " << shell->getSize() << endl; //<< " " << FortranFormat(shell->getScale()) << endl;
                                    int _sh_idx =0;
                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                        _sh_idx++;
                                        _el_file << " " << _sh_idx << " " << FortranFormat(gaussian->decay);
                                        for (unsigned _icontr = 0; _icontr < gaussian->contraction.size(); _icontr++) {
                                            if (gaussian->contraction[_icontr] != 0.0) {
                                                _el_file << " " << FortranFormat(gaussian->contraction[_icontr]);
                                            }
                                        }
                                        _el_file << endl;
                                    }
                                }
                            }
                        }
                    }
                    _el_file << "STOP\n";
                    _el_file.close();
                
                
                    _com_file << "%basis\n "  << endl;
                    _com_file << "GTOName" << " " << "=" << "\"system.bas\";"  << endl;
                } // write_basis set
                
                /* WRITING ECP INTO system.inp FILE for ORCA**/
                if ( _write_pseudopotentials ){
                    std::string pseudopotential_name("ecp");
                    _com_file << endl;
                    Angularmom lmaxnum; //defining lmaxnum for lmaxNum2lmaxName
                    list<std::string> elements;
                    elements.push_back("H");
                    elements.push_back("He");
                    BasisSet ecp;
                    ecp.LoadPseudopotentialSet(pseudopotential_name);
                    LOG(CTP::logDEBUG, *_pLog) << "Loaded Pseudopotentials " << pseudopotential_name << flush;
                    for (sit = segments.begin(); sit != segments.end(); ++sit) {
                        std::vector< CTP::Atom* > atoms = (*sit)-> Atoms();
                        std::vector< CTP::Atom* >::iterator it;
                        for (it = atoms.begin(); it < atoms.end(); it++) {
                            std::string element_name = (*it)->getElement();
                            list<std::string>::iterator ite;
                            ite = find(elements.begin(), elements.end(), element_name);
                            if (ite == elements.end()) {
                                elements.push_back(element_name);
                                Element* element = ecp.getElement(element_name);
                                _com_file << "\n" << "NewECP" << " " <<  element_name  << endl;
                                _com_file << "N_core" << " " << element->getNcore()  << endl;
                                //lmaxnum2lmaxname 
                                _com_file << "lmax" << " " << lmaxnum.getLName(element->getLmax()) << endl;
                             
                                //For Orca the order doesn't matter but let's write it in ascending order
                                // write remaining shells in ascending order s,p,d...
                                for (int i = 0; i < element->getLmax(); i++) {
                                    for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                        Shell* shell = (*its);
                                        if (shell->getLmax() == i) {
                                        // shell type, number primitives, scale factor
                                        _com_file << shell->getType() << " " <<  shell->getSize() << endl;
                                        int _sh_idx=0;
                                        for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                            GaussianPrimitive* gaussian = *itg;
                                            _sh_idx++;
                                            _com_file << _sh_idx << " " << gaussian->decay << " " << gaussian->contraction[0] << " " << gaussian->power << endl;
                                        }
                                    }
                                }
                            }
                                
                            for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                Shell* shell = (*its);
                                // shell type, number primitives, scale factor
                                if (shell->getLmax() == element->getLmax()) {
                                    _com_file << shell->getType() << " " <<  shell->getSize() << endl;
                                    // _com_file << shell->getSize() << endl;
                                    int _sh_idx=0;
                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                             _sh_idx++;
                                        _com_file << _sh_idx << " " << gaussian->decay << " " << gaussian->contraction[0] << " " << gaussian->power << endl;
                                        }
                                    }
                                }
                                
                               _com_file << "end\n "  <<  "\n" << endl; 
                            }
                           
                        }
                    }  
                        
        }//if(_write_pseudopotentials)     
        /* END   OF WRITING ECP INTO system.inp FILE for ORCA*************/                
        _com_file << "end\n "  <<  "\n" << endl;   //This end is for the basis set block
                
    } else {
            std::vector< CTP::QMAtom* > *qmatoms = orbitals_guess->getAtoms();
            std::vector< CTP::QMAtom* >::iterator it;

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
            _com_file << "* \n" << endl;
            _com_file << "%pal\n "  <<  "nprocs " <<  _threads  << "\nend" << "\n" << endl;

            
            // basis set info
            if ( _write_basis_set) {
                Elements _elements;
   
                list<std::string> elements;
                BasisSet bs;
                bs.LoadBasisSet(_basisset_name);
                LOG(CTP::logDEBUG, *_pLog) << "Loaded Basis Set " << _basisset_name << flush;
                ofstream _el_file;
                std::string _el_file_name = _run_dir + "/" +  "system.bas";
                _el_file.open(_el_file_name.c_str());
                //_com_file << "@" << "system.bas" << endl;
                _el_file << "$DATA" << endl;
                                                    
                   for (it = qmatoms->begin(); it < qmatoms->end(); it++) {
                        if (!(*it)->from_environment) {
                            std::string element_name = (*it)->type;
                            list<std::string>::iterator ite;
                            ite = find(elements.begin(), elements.end(), element_name);
                            if (ite == elements.end()) {
                                elements.push_back(element_name);
                                Element* element = bs.getElement(element_name);
                                _el_file << _elements.getEleFull(element_name) << endl; 
                                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                    Shell* shell = (*its);
                                    _el_file  << shell->getType() << " " << shell->getSize() << endl; //<< " " << FortranFormat(shell->getScale()) << endl;
                                    int _sh_idx =0;
                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                        _sh_idx++;
                                        _el_file << " " << _sh_idx << " " << FortranFormat(gaussian->decay);
                                        for (unsigned _icontr = 0; _icontr < gaussian->contraction.size(); _icontr++) {
                                            if (gaussian->contraction[_icontr] != 0.0) {
                                                _el_file << " " << FortranFormat(gaussian->contraction[_icontr]);
                                            }
                                        }
                                        _el_file << endl;
                                    }
                                }
                            }
                        }
                    }
                    _el_file << "STOP\n";
                    _el_file.close();
                
                
                    _com_file << "%basis\n "  << endl;
                    _com_file << "GTOName" << " " << "=" << "\"system.bas\";"  << endl;
                } // write_basis set
            
            // ECPs
               /* WRITING ECP INTO system.inp FILE for ORCA**/
                if ( _write_pseudopotentials ){
                    std::string pseudopotential_name("ecp");
                    _com_file << endl;
                    Angularmom lmaxnum; //defining lmaxnum for lmaxNum2lmaxName
                    list<std::string> elements;
                    elements.push_back("H");
                    elements.push_back("He");
                    BasisSet ecp;
                    ecp.LoadPseudopotentialSet(pseudopotential_name);
                    LOG(CTP::logDEBUG, *_pLog) << "Loaded Pseudopotentials " << pseudopotential_name << flush;
                    
                    
                    
                     for (it = qmatoms->begin(); it < qmatoms->end(); it++) {
                        if (!(*it)->from_environment) {
                            std::string element_name = (*it)->type;
                            
                            
                            list<std::string>::iterator ite;
                            ite = find(elements.begin(), elements.end(), element_name);
                            if (ite == elements.end()) {
                                elements.push_back(element_name);
                                Element* element = ecp.getElement(element_name);
                                _com_file << "\n" << "NewECP" << " " <<  element_name  << endl;
                                _com_file << "N_core" << " " << element->getNcore()  << endl;
                                //lmaxnum2lmaxname 
                                _com_file << "lmax" << " " << lmaxnum.getLName(element->getLmax()) << endl;
                             
                                //For Orca the order doesn't matter but let's write it in ascending order
                                // write remaining shells in ascending order s,p,d...
                                for (int i = 0; i < element->getLmax(); i++) {
                                    for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                        Shell* shell = (*its);
                                        if (shell->getLmax() == i) {
                                        // shell type, number primitives, scale factor
                                        _com_file << shell->getType() << " " <<  shell->getSize() << endl;
                                        int _sh_idx = 0;
                                        for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                            GaussianPrimitive* gaussian = *itg;
                                            _sh_idx++;
                                            _com_file << _sh_idx << " " << gaussian->decay << " " << gaussian->contraction[0] << " " << gaussian->power << endl;
                                        }
                                    }
                                }
                            }
                                
                            for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                                Shell* shell = (*its);
                                // shell type, number primitives, scale factor
                                if (shell->getLmax() == element->getLmax()) {
                                    _com_file << shell->getType() << " " <<  shell->getSize() << endl;
                                    // _com_file << shell->getSize() << endl;
                                    int _sh_idx = 0;
                                    for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                                        GaussianPrimitive* gaussian = *itg;
                                        _sh_idx++;
                                        _com_file << _sh_idx << " " << gaussian->decay << " " << gaussian->contraction[0] << " " << gaussian->power << endl;
                                        }
                                    }
                                }
                                
                               _com_file << "end\n "  <<  "\n" << endl; 
                            }
                           
                        }
                    }
                } // write pseudopotentials
                
                
                /* END   OF WRITING ECP INTO system.inp FILE for ORCA*************/                
                _com_file << "end\n "  <<  "\n" << endl;   //This end is for the basis set block
             
            
            
            // actual writing of charges    
            _crg_file.open( _crg_file_name_full.c_str() ); 
            int _total_background = 0;
            for (it = qmatoms->begin(); it < qmatoms->end(); it++) {
                if ((*it)->from_environment) {
                    if ((*it)->charge != 0.0) _total_background++;
                }
            }

            _crg_file << _total_background << endl;
            for (it = qmatoms->begin(); it < qmatoms->end(); it++) {
                if ((*it)->from_environment) {
                    boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
                    fmt % (*it)->charge % (*it)->x % (*it)->y % (*it)->z ;
                    if ((*it)->charge != 0.0) _crg_file << fmt << endl;
                }
            }
                
         } // if (! _write_charges )
    

    
    
 
   
   /*************************************WRITING BASIS SET INTO system.bas FILE**********************************/
 
    /***************************** END OF WRITING BASIS SET INTO system.bas FILE**********************************************/
 
     
    // writing scratch_dir info
    if ( _scratch_dir != "" ) {
        
        LOG(CTP::logDEBUG,*_pLog) << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

        // boost::filesystem::create_directories( _scratch_dir + temp_suffix );
        std::string _temp( "scratch_dir " + _scratch_dir + temp_suffix + "\n" );
        //_com_file << _temp ;
    }

    _com_file << _options << "\n";

    
    if ( _write_guess ) { 
       
            throw std::runtime_error( "Not implemented in orca");
    }   

    _com_file << endl;
    _com_file.close();
    
    
    
      // and now generate a shell script to run both jobs, if neccessary
    LOG(CTP::logDEBUG, *_pLog) << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

    _scratch_dir = scratch_dir_backup + temp_suffix;
            
            //boost::filesystem::create_directories(_scratch_dir + temp_suffix);
            //std::string _temp("scratch_dir " + _scratch_dir + temp_suffix + "\n");
            //_com_file << _temp;
    WriteShellScript();
    _scratch_dir = scratch_dir_backup;
   
    return true;
}

bool Orca::WriteShellScript() {
    ofstream _shell_file;
    
    std::string _shell_file_name_full = _run_dir + "/" + _shell_file_name;
            
    _shell_file.open ( _shell_file_name_full.c_str() );

    _shell_file << "#!/bin/bash" << endl ;
    _shell_file << "mkdir -p " << _scratch_dir << endl;
    
    if ( _threads == 1 ){ 
        _shell_file << _executable << " " << _input_file_name << " > " <<  _log_file_name  << endl;//" 2> run.error" << endl;    
    } else {
        _shell_file << _executable << " "<< _input_file_name << " > " <<  _log_file_name <<endl;// " 2> run.error" << endl;    
    }
    _shell_file.close();
    
    return true;   
}

/**
 * Runs the Orca job. 
 */
bool Orca::Run()
{

    LOG(CTP::logDEBUG,*_pLog) << "Running Orca job" << flush;
    
    if (std::system(NULL)) {
        
        // Orca overrides input information, if *.db and *.movecs files are present
        // better trash the old version
        std::string file_name = _run_dir + "/system.db";
        remove ( file_name.c_str() );
        file_name = _run_dir + "/" + _log_file_name;
        remove ( file_name.c_str() );
        file_name = _run_dir + "/" + _orb_file_name;
        //remove ( file_name.c_str() );
               
        std::string _command;
        if ( _threads == 1 ) {
            _command = "cd " + _run_dir + "; sh " + _shell_file_name;
        } else {
            _command = "cd " + _run_dir + "; sh " + _shell_file_name;
        }
        //LOG(CTP::logDEBUG,*_pLog) << _command << flush;
       int check=std::system(_command.c_str());
        if (check==-1){
            LOG(CTP::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
            return false;
        }
       
        if ( CheckLogFile() ) {
            LOG(CTP::logDEBUG,*_pLog) << "Finished Orca job" << flush;
            return true;
        } else {
            LOG(CTP::logDEBUG,*_pLog) << "Orca job failed" << flush;
        }
    }
    else {
        LOG(CTP::logERROR,*_pLog) << _input_file_name << " failed to start" << flush; 
        return false;
    }
    
    return true;
}

/**
 * Cleans up after the Orca job
 */
void Orca::CleanUp() {
    
    // cleaning up the generated files
    if ( _cleanup.size() != 0 ) {
        Tokenizer tok_cleanup(_cleanup, ",");
        std::vector <std::string> _cleanup_info;
        tok_cleanup.ToVector(_cleanup_info);
        
        std::vector<std::string> ::iterator it;
               
        for (it = _cleanup_info.begin(); it != _cleanup_info.end(); ++it) {
            if ( *it == "inp" ) {
                std::string file_name = _run_dir + "/" + _input_file_name;
                remove ( file_name.c_str() );
            }
            
            if ( *it == "bas" ) {
                std::string file_name = _run_dir +"/system.bas";
                remove ( file_name.c_str() );
            }
            
            if ( *it == "log" ) {
                std::string file_name = _run_dir + "/" + _log_file_name;
                remove ( file_name.c_str() );
            }

           if ( *it == "gbw" ) {
                std::string file_name = _run_dir + "/" + _orb_file_name;
                remove ( file_name.c_str() );
            }
            
            if ( *it == "ges" ) {
                std::string file_name = _run_dir + "/system.ges" ;
                remove ( file_name.c_str() );
            }    
            if ( *it == "prop" ) {
                std::string file_name = _run_dir + "/system.prop" ;
                remove ( file_name.c_str() );
            }     
        }
    }
    
}



/**
 * Reads in the MO coefficients from an Orca gbw file
 */
bool Orca::ParseLogFile( Orbitals* _orbitals )
{
    static const double _conv_Hrt_eV = 27.21138386;
    
    _orbitals->setQMpackage("orca");
    
    LOG(CTP::logDEBUG,*_pLog) << "Parsing " << _log_file_name << flush;
      // return true;
  std::string _log_file_name_full =  _run_dir + "/" + _log_file_name;
   // check if LOG file is complete
   if ( !CheckLogFile() ) return false;
    
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
    

    
    std::ifstream _input_file( _log_file_name_full.c_str() );
    
    if (_input_file.fail()) {
        LOG( CTP::logERROR, *_pLog ) << "File " << _log_file_name_full << " not found " << flush;
        return false;
    } else {
        LOG(CTP::logDEBUG, *_pLog) << "Reading Coordinates and occupationnumbers and energies from " << _log_file_name_full << flush;
    }
    
      
          //Coordinates of the final configuration depending on whether it is an optimization or not
         

    
        
        
    while (_input_file) {
        getline(_input_file,_line);
        boost::trim(_line);
        
        
        
        if ( _is_optimization ){
              throw runtime_error("Not implemented yet!");
        }
        bool _found_optimization=true;
        
        std::string::size_type coordinates_pos = _line.find("CARTESIAN COORDINATES (ANGSTROEM)");
        
        if ( _found_optimization && coordinates_pos != std::string::npos) {
            LOG(CTP::logDEBUG,*_pLog) << "Getting the coordinates" << flush;
            
            //_has_coordinates = true;
            bool _has_QMAtoms = _orbitals->hasQMAtoms();

            // three garbage lines
            getline(_input_file, _line);
            // now starts the data in format
            // _id type Qnuc x y z 
            std::vector<std::string> _row;
            getline(_input_file, _line);
            boost::trim( _line );

            boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on); 
            int nfields =  _row.size();

            int atom_id = 0;   
            while ( nfields == 4 ) {
                //int atom_id = boost::lexical_cast< int >( _row.at(0) );
                //int atom_number = boost::lexical_cast< int >( _row.at(0) );
                std::string _atom_type = _row.at(0);
                double _x =  boost::lexical_cast<double>( _row.at(1) );
                double _y =  boost::lexical_cast<double>( _row.at(2) );
                double _z =  boost::lexical_cast<double>( _row.at(3) );
                //if ( tools::globals::verbose ) cout << "... ... " << atom_id << " " << atom_type << " " << atom_charge << endl;
                getline(_input_file, _line);
                boost::trim( _line );
                boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);  
                nfields =  _row.size();

                if ( _has_QMAtoms == false ) {
                    _orbitals->AddAtom( _atom_type, _x, _y, _z );
                } else {
                                       
                    CTP::QMAtom* pAtom = _orbitals->_atoms.at( atom_id  );
                    pAtom->type = _atom_type;
                    pAtom->x = _x;
                    pAtom->y = _y;
                    pAtom->z = _z;
                    atom_id++;
                }
                    
            }
      
           

        }
        
        
        std::string::size_type energy_pos = _line.find("Total Energy");
       if (energy_pos != std::string::npos) {
            //cout << _line << endl;
            boost::algorithm::split(results, _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
            std::string _energy = results[3];
            boost::trim( _energy );
            //cout << _energy << endl; 
            _orbitals->setQMEnergy ( _conv_Hrt_eV * boost::lexical_cast<double>(_energy) );
            LOG(CTP::logDEBUG, *_pLog) << (boost::format("QM energy[eV]: %4.6f ") % _orbitals->getQMEnergy()).str() << flush;
            // _orbitals->_has_qm_energy = true;
            }
        
        /* Check for ScaHFX = factor of HF exchange included in functional */
              std::string::size_type HFX_pos = _line.find("Fraction HF Exchange ScalHFX");
              if (HFX_pos != std::string::npos) {
                  boost::algorithm::split(results, _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
                  double _ScaHFX = boost::lexical_cast<double>(results.back());
                  _orbitals->setScaHFX(_ScaHFX);
                  LOG(CTP::logDEBUG, *_pLog) << "DFT with " << _ScaHFX << " of HF exchange!" << flush;
              }
     
    //Finding Basis Dimension, the number of energy levels
        std::string::size_type dim_pos = _line.find("Basis Dimension");
        if (dim_pos != std::string::npos) {
                
            boost::algorithm::split(results, _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
            //_has_basis_dim = true;
            std::string _dim = results[4];  //The 4th element of results vector is the Basis Dim
            boost::trim( _dim );
            _levels = boost::lexical_cast<int>(_dim);
            //cout <<  boost::lexical_cast<int>(_dim) << endl;
            //_basis_size = _levels;
              LOG(CTP::logDEBUG,*_pLog) << "Basis Dimension: " << _levels << flush;
              LOG( CTP::logDEBUG, *_pLog ) << "Energy levels: " << _levels << flush;
        }
        /********************************************************/
        
        
        
        
        
        
        std::string::size_type OE_pos = _line.find("ORBITAL ENERGIES");
        if (OE_pos != std::string::npos) {
        getline(_input_file,_line);
        getline(_input_file,_line);
        getline(_input_file,_line);
        if (_line.find("E(Eh)")==std::string::npos){
            LOG(CTP::logDEBUG,*_pLog) << "Warning: Orbital Energies not found in log file" << flush;
        }
        for (unsigned i=0; i<_levels ; i++){
            getline(_input_file,_line);
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
           
                    
            std::string _no =results[0];
            
            boost::trim( _no );
            unsigned levelnumber= boost::lexical_cast<unsigned>(_no);
            if (levelnumber!=i){
                LOG(CTP::logDEBUG,*_pLog) << "Have a look at the orbital energies something weird is going on" << flush;
            }
            std::string _oc = results[1];
            boost::trim( _oc );
            double occ =boost::lexical_cast<double>(_oc);
            // We only count alpha electrons, each orbital must be empty or doubly occupied
            if (occ==2){
                _number_of_electrons++;
                _occ[i]=occ;
            }
            else if (occ==0) {
                _occ[i]=occ;
            }
            else {
                throw runtime_error("Only empty or doubly occupied orbitals are allowed not running the right kind of DFT calculation");
            }
            
            std::string _e = results[2];
            boost::trim( _e );
            _energies [i] = boost::lexical_cast<double>(_e) ;
        }
        }
            /*
                 *  Partial charges from the input file
                 */
                std::string::size_type charge_pos = _line.find("CHELPG Charges");

                if (charge_pos != std::string::npos && _get_charges) {
                    LOG(CTP::logDEBUG, *_pLog) << "Getting charges" << flush;
                    //_has_charges = true;
                    getline(_input_file, _line);
                    //getline(_input_file, _line);

                    bool _has_atoms = _orbitals->hasQMAtoms();

                    std::vector<std::string> _row;
                    getline(_input_file, _line);
                    boost::trim(_line);
                    //cout << _line << endl;
                    boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    int nfields = _row.size();
                    //cout << _row.size() << endl;

                    while (nfields == 4) {
                        int atom_id = boost::lexical_cast< int >(_row.at(0));
                        atom_id++;
                        //int atom_number = boost::lexical_cast< int >(_row.at(0));
                        std::string atom_type = _row.at(1);
                        double atom_charge = boost::lexical_cast< double >(_row.at(3));
                        if ( votca::tools::globals::verbose ) cout << "... ... " << atom_id << " " << atom_type << " " << atom_charge << endl;
                        getline(_input_file, _line);
                        boost::trim(_line);
                        boost::algorithm::split(_row, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                        nfields = _row.size();

                        if (_has_atoms == false) {
                            _orbitals->AddAtom(atom_type, 0, 0, 0, atom_charge);
                        } else {
                            CTP::QMAtom* pAtom = _orbitals->_atoms.at(atom_id - 1);
                            pAtom->type = atom_type;
                            pAtom->charge = atom_charge;
                        }

                    }
                    //_orbitals->_has_atoms = true;
                }
        }

    
    LOG(CTP::logDEBUG,*_pLog) << "Alpha electrons: " << _number_of_electrons << flush ;
    int _occupied_levels = _number_of_electrons;
    int _unoccupied_levels = _levels - _occupied_levels;
    LOG(CTP::logDEBUG,*_pLog) << "Occupied levels: " << _occupied_levels << flush;
    LOG(CTP::logDEBUG,*_pLog) << "Unoccupied levels: " << _unoccupied_levels << flush;  
   
        
   /************************************************************/
 
    // copying information to the orbitals object
   /* _orbitals->setBasisSetSize(  _basis_size );*/
    _orbitals->setBasisSetSize(  _levels );
   
    _orbitals->setNumberOfElectrons( _number_of_electrons );

    _orbitals->setNumberOfLevels( _occupied_levels , _unoccupied_levels );
    
    _orbitals->setSelfEnergy(0.0);
    
   // copying energies to a matrix  
   _orbitals->_mo_energies.resize( _levels );
   //_level = 1;
   for(size_t i=0; i < _orbitals->_mo_energies.size(); i++) {
         _orbitals->_mo_energies[i] = _energies[ i ];
    }

   
   // copying orbitals to the matrix
   
 
   
   
   //cout << _mo_energies << endl;   
   // cout << _mo_coefficients << endl; 
   
   // cleanup
  // _coefficients.clear();
   _energies.clear();
   _occ.clear();
   
   
   
   

   
   
   
   LOG(CTP::logDEBUG, *_pLog) << "Done reading Log file" << flush;

   return true;
}//ParseOrbitalFile(Orbital* _orbital)

bool Orca::CheckLogFile() {
    
    // check if the log file exists
    ifstream _input_file((_run_dir + "/" + _log_file_name).c_str());
    
    if (_input_file.fail()) {
        LOG(CTP::logERROR,*_pLog) << "Orca LOG is not found" << flush;
        return false;
    };
    
    std::string _line;
    while (_input_file) {
        getline(_input_file,_line);
        boost::trim(_line);
        

        std::string::size_type error = _line.find("FATAL ERROR ENCOUNTERED");
        
        if (error != std::string::npos) {
            LOG(CTP::logERROR,*_pLog) << "ORCA encountered a fatal error, maybe a look in the log file may help." << flush;
            return false;
        } 
        error = _line.find("mpirun detected that one or more processes exited with non-zero status");
        
        if (error != std::string::npos) {
            LOG(CTP::logERROR,*_pLog) << "ORCA had an mpi problem, maybe your openmpi version is not good." << flush;
            return false;
        } 
    } 
    return true;
}

 // Parses the Orca gbw file and stores data in the Orbitals object 
 
bool Orca::ParseOrbitalsFile( Orbitals* _orbitals ) {
if ( !CheckLogFile() ) return false;
std::vector<double> _coefficients;
int _basis_size=_orbitals->getBasisSetSize();
int _levels=_orbitals->getNumberOfLevels();

if (_basis_size==0 || _levels==0){
    throw runtime_error("Basis size not set, calculator does not parse log file first");
}

   LOG(CTP::logDEBUG,*_pLog) << "Reading the gbw file, this may or may not work so be careful: " << flush ;
    ifstream infile;
    infile.open(( _run_dir + "/" + _orb_file_name).c_str(), ios::binary | ios::in);
    if (!infile){
        throw runtime_error("Could not open "+_orb_file_name+" file");
    }
    infile.seekg(24, ios::beg); 
    char* buffer= new char [8];
    infile.read(buffer,8);
    long int offset= *((long int*)buffer);
    
    infile.seekg(offset,ios::beg);
    infile.read(buffer,4);
    int op_read= *((int*)buffer);
    infile.seekg(offset+4,ios::beg);
    infile.read(buffer,4);
    int dim_read= *((int*)buffer);
    infile.seekg(offset+8,ios::beg);
    LOG(CTP::logDEBUG,*_pLog) << "Number of operators: "<< op_read<< " Basis dimension: "<<dim_read<< flush ;
    int n=op_read*dim_read*dim_read;
    delete[] buffer;
    buffer =new char [8];
    for (int i=0;i<n;i++){
        infile.read(buffer,8);
        double mocoeff=*((double*)buffer);
         //LOG(CTP::logDEBUG,*_pLog) << mocoeff<< flush ;
        _coefficients.push_back(mocoeff);
    }
    delete[] buffer;
    
    infile.close();
    //cout<< "basissize " <<_basis_size << endl;
    //cout << "coeffvektor size "<< _coefficients.size() << endl;
    
    
    // i -> MO, j -> AO 
      (_orbitals->_mo_coefficients).resize( _levels, _basis_size );     
   for(size_t i = 0; i < _orbitals->_mo_coefficients.size1(); i++) {
      for(size_t j = 0 ; j < _orbitals->_mo_coefficients.size2(); j++) {
         _orbitals->_mo_coefficients(i,j) = _coefficients[j*_basis_size+i];
         //cout <<  _coefficients[i][j] << endl;
         //cout << i << " " << j << endl;
        }
    }
   
      //cout<<"MO1:1 :" <<setprecision(15)<<_orbitals->_mo_coefficients(0,0)<< endl;
      //cout<<"MO-1:-1"<<setprecision(15)<< _orbitals->_mo_coefficients(_levels-1,_levels-1)<< endl;
      
      
    LOG(CTP::logDEBUG,*_pLog) << "Done parsing" << flush;
    return true;
}



/**
 * Converts the Orca data stored in the Orbitals object to GW input format
 */
bool Orca::ConvertToGW( Orbitals* _orbitals ) {
    cerr << "Tried to convert to GW from Orca package. ";
    throw std::runtime_error( "Conversion not implemented yet!");
}


std::string Orca::FortranFormat( const double &number ) {
    std::stringstream _ssnumber;
    if ( number >= 0) {
        _ssnumber << "    ";
    } else{
        _ssnumber << "   ";
    }
        
    _ssnumber <<  setiosflags(ios::fixed) << setprecision(15) << std::scientific << number;
    std::string _snumber = _ssnumber.str(); 
    //boost::replace_first(_snumber, "e", "D");
    return _snumber;
}
        



}}
