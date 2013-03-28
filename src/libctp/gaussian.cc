/*
 *            Copyright 2009-2012 The VOTCA Development Team
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

#include "votca/ctp/gaussian.h"
#include "votca/ctp/segment.h"
#include <stdio.h>
#include <iomanip>
#include <sys/stat.h>

using namespace std;

namespace votca { namespace ctp {

Gaussian::Gaussian( tools::Property *opt ) { 
    
    string key = "package";

    string _name = opt->get(key+".name").as<string> ();
    
    if ( _name != "gaussian" ) {
        cerr << "Tried to use " << _name << " package. ";
        throw std::runtime_error( "Package is not supported.");
    }
    
    _executable =       opt->get(key + ".executable").as<string> ();
    _functional =       opt->get(key + ".functional").as<string> ();
    _basis_set =        opt->get(key + ".basisset").as<string> ();
    _charge =           opt->get(key + ".charge").as<int> ();
    _spin =             opt->get(key + ".spin").as<int> ();
    _options =          opt->get(key + ".options").as<string> ();
    _memory =           opt->get(key + ".memory").as<string> ();
    _threads =          opt->get(key + ".threads").as<int> ();
    _chk_file_name  =   opt->get(key + ".checkpoint").as<string> ();
    _scratch_dir =      opt->get(key + ".scratch").as<string> ();
    _cleanup =          opt->get(key + ".cleanup").as<string> ();
        
};   
    
Gaussian::~Gaussian() { 
}  

/**
 * Prepares the com file from a vector of segments
 */
bool Gaussian::WriteInputFile( Segment *seg ) {

    vector< Atom* > ::iterator ait;

    int qmatoms = 0;

    vector< Atom* > _atoms;
    _atoms  = seg-> Atoms();

    ofstream _com_file;
    
    string _com_file_name_full = _run_dir + "/" + _com_file_name;
    
    _com_file.open ( _com_file_name_full.c_str() );
    // header 
    if ( _chk_file_name.size() != 0 ) {
        _com_file << "%chk = " << _chk_file_name << endl;
    }
    
    if ( _memory.size() != 0 ) {
        _com_file << "%mem = " << _memory << endl ;
    }
    
    if ( _threads != 0 ) {
         _com_file << "%nprocshared = "  << _threads << endl;
    }
    
    if ( _functional.size() != 0 and _basis_set.size() != 0 ) {
        _com_file <<  "# " <<  _functional << "/" <<  _basis_set;
    }
    
    if ( _options.size() != 0 ) {
        _com_file <<  "  " << _options << endl ;
    }
    
    _com_file << endl << seg->getName() << endl << endl;
    _com_file << setw(2) << _charge << setw(2) << _spin << endl;
    
    for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

        if ((*ait)->HasQMPart() == false) { continue; }

        vec     pos = (*ait)->getQMPos();
        string  name = (*ait)->getElement();

        //fprintf(out, "%2s %4.7f %4.7f %4.7f \n"
        _com_file << setw(3) << name.c_str() 
                  << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX()*10
                  << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY()*10
                  << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ()*10 
                  << endl;
    }
    
    _com_file << endl;
    _com_file.close();
}

bool Gaussian::WriteShellScript() {
    ofstream _shell_file;
    
    string _shell_file_name_full = _run_dir + "/" + _shell_file_name;
            
    _shell_file.open ( _shell_file_name_full.c_str() );

    _shell_file << "#!/bin/tcsh" << endl ;
    _shell_file << "mkdir -p " << _scratch_dir << endl;
    _shell_file << "setenv GAUSS_SCRDIR " << _scratch_dir << endl;
    _shell_file << _executable << " " << _com_file_name << endl;    
    _shell_file.close();
}

/**
 * Runs the Gaussian job
 */
bool Gaussian::Run()
{

    if (system(NULL)) {
        // if scratch is provided, run the shell script; 
        // otherwise run gaussian directly and rely on global variables 
        string _command;
        if ( _scratch_dir.size() != 0 ) {
            _command  = "cd " + _run_dir + "; tcsh " + _shell_file_name;
        }
        else {
            _command  = "cd " + _run_dir + "; " + _executable + " " + _com_file_name;
        }
        
        int i = system ( _command.c_str() );
    }
    else {
        cerr << "The job " << _com_file_name << " failed to complete" << endl; 
        exit (EXIT_FAILURE);
    }

}

/**
 * Cleans up after the Gaussian job
 */
void Gaussian::CleanUp( string ID ) {
    
    // cleaning up the generated files
    if ( _cleanup.size() != 0 ) {
        Tokenizer tok_cleanup(_cleanup, " \t\n");
        vector <string> _cleanup_info;
        tok_cleanup.ToVector(_cleanup_info);
        
        vector<string> ::iterator it;
        
        for (it = _cleanup_info.begin(); it != _cleanup_info.end(); ++it) {
            if ( *it == "xyz" || *it == "com" || *it == "log" ) { 
                string file_name = _run_dir + "/mol_" + ID + "." + *it;
                remove ( file_name.c_str() );
            }
            if ( *it == "fort.7" ) {
                string file_name = _run_dir + "/" + *it;
                remove ( file_name.c_str() );
            }            
        }
    }
    
}


bool Gaussian::ParseLog( Orbitals* _orbitals ) {
    
    string _line;
    unsigned _basis_size = 0;
    bool _read_overlap = false;
    bool _verbose = true;
    
    if ( _verbose )  cout << endl << "... ... Parsing " << _com_file_name << endl;
    
    ifstream _input_file( _com_file_name.c_str() );
    if (_input_file.fail()) {
        throw std::runtime_error("COM file is not found." );
        return 1;
    };

    while ( _input_file  ) {

        getline(_input_file, _line);
        
        // number of electrons
        std::string::size_type electrons_pos = _line.find("alpha electrons");
        if ( electrons_pos != std::string::npos ) cout << _input_file.tellg();
               
        // basis set size
        std::string::size_type basis_pos = _line.find("basis functions");
        if ( basis_pos != std::string::npos ) cout << _input_file.tellg();
        
        // energies of occupied/unoccupied levels
        std::string::size_type eigenvalues_pos = _line.find("Alpha");
         if ( eigenvalues_pos != std::string::npos ) cout << _input_file.tellg();
       
    }    
}


}}
