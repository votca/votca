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

#include "gw.h"
#include "votca/ctp/segment.h"
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

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
void GW::Initialize( Property *options ) {

     // NWChem file names
    string fileName = "system";

    _xyz_file_name = fileName + ".xyz";
    _input_file_name = fileName + ".nw";
    _log_file_name = fileName + ".log"; 
    _orb_file_name = fileName + ".movecs" ;               

    string key = "package";
    string _name = options->get(key+".name").as<string> ();
    
    if ( _name != "gw" ) {
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
    
}    

/**
 * Prepares the *.nw file from a vector of segments
 * Appends a guess constructed from monomer orbitals if supplied
 */
bool GW::WriteInputFile( vector<Segment* > segments, Orbitals* orbitals_guess )
{
         
}

/**
 * Runs the NWChem job. 
 */
bool GW::Run()
{

    LOG(logDEBUG,*_pLog) << "Running NWChem job" << flush;
    
    if (system(NULL)) {
        
        // NWChem overrides input information, if *.db and *.movecs files are present
        // better trash the old version
        string file_name = _run_dir + "/system.db";
        remove ( file_name.c_str() );
        file_name = _run_dir + "/" + _log_file_name;
        remove ( file_name.c_str() );
        file_name = _run_dir + "/" + _orb_file_name;
        //remove ( file_name.c_str() );
        
                
        // if threads is provided and > 1, run mpi; 
        string _command;
        if ( _threads == 1 ) {
            _command  = "cd " + _run_dir + "; " + _executable + " " + _input_file_name + "> " +  _log_file_name ;
        } else {
            _command  = "cd " + _run_dir + "; mpirun -np " +  boost::lexical_cast<string>(_threads) + " " + _executable + " " + _input_file_name + "> "+  _log_file_name ;
        }
        
        int i = system ( _command.c_str() );
        
        if ( CheckLogFile() ) {
            LOG(logDEBUG,*_pLog) << "Finished NWChem job" << flush;
            return true;
        } else {
            LOG(logDEBUG,*_pLog) << "NWChem job failed" << flush;
        }
    }
    else {
        LOG(logERROR,*_pLog) << _input_file_name << " failed to start" << flush; 
        return false;
    }
    



}

/**
 * Cleans up after the NWChem job
 */
void GW::CleanUp() {
    
    // cleaning up the generated files
    if ( _cleanup.size() != 0 ) {
        Tokenizer tok_cleanup(_cleanup, " \t\n");
        vector <string> _cleanup_info;
        tok_cleanup.ToVector(_cleanup_info);
        
        vector<string> ::iterator it;
               
        for (it = _cleanup_info.begin(); it != _cleanup_info.end(); ++it) {
            if ( *it == "nw" ) {
                string file_name = _run_dir + "/" + _input_file_name;
                remove ( file_name.c_str() );
            }
            
            if ( *it == "db" ) {
                string file_name = _run_dir + "/system.db";
                remove ( file_name.c_str() );
            }
            
            if ( *it == "log" ) {
                string file_name = _run_dir + "/" + _log_file_name;
                remove ( file_name.c_str() );
            }

           if ( *it == "movecs" ) {
                string file_name = _run_dir + "/system.movecs";
                remove ( file_name.c_str() );
            }
            
            if ( *it == "gridpts" ) {
                string file_name = _run_dir + "/system.gridpts.*" ;
                remove ( file_name.c_str() );
            }            
        }
    }
    
}



/**
 * Reads in the MO coefficients from a NWChem movecs file
 */
bool GW::ParseOrbitalsFile( Orbitals* _orbitals )
{
   return true;
}

bool GW::CheckLogFile() {
    
    // check if the log file exists
    char ch;
    ifstream _input_file((_run_dir + "/" + _log_file_name).c_str());
    
    if (_input_file.fail()) {
        LOG(logERROR,*_pLog) << "NWChem LOG is not found" << flush;
        return false;
    };

    
    /* Checking the log file is a pain in the *** since NWChem throws an error
     * for our 'iterations 1'  runs (even though it outputs the required data
     * correctly. The only way that works for both scf and noscf runs is to
     * check for "Total DFT energy" near the end of the log file. 
     */
    
    _input_file.seekg(0,ios_base::end);   // go to the EOF

    std::string::size_type total_energy_pos = std::string::npos;
    std::string::size_type diis_pos = std::string::npos;
    do {
       // get the beginning of the line 
       do {
          _input_file.seekg(-2,ios_base::cur);
          _input_file.get(ch);   
          //cout << "\nNext Char: " << ch << " TELL G " <<  (int)_input_file.tellg() << endl;
       } while ( ch != '\n' || (int)_input_file.tellg() == -1 );
            
       string _line;            
       getline(_input_file,_line);                      // Read the current line
       // cout << "\nResult: " << _line << '\n';     // Display it
       total_energy_pos = _line.find("Total DFT energy");
       diis_pos = _line.find("diis");
       // whatever is found first, determines the completeness of the file
       if (total_energy_pos != std::string::npos) {
          return true;
       } else if (diis_pos != std::string::npos) {
           LOG(logERROR,*_pLog) << "NWChem LOG is incomplete" << flush;
           return false;
       } else {
           // go to previous line
           //_input_file.get(ch); 
           do {
           _input_file.seekg(-2,ios_base::cur);
           _input_file.get(ch);   
           //cout << "\nChar: " << ch << endl;
         } while ( ch != '\n' || (int)_input_file.tellg() == -1);
       }
    } while ( total_energy_pos == std::string::npos && diis_pos == std::string::npos );
    
    
    _input_file.close();
    
}

/**
 * Parses the Gaussian Log file and stores data in the Orbitals object 
 */
bool GW::ParseLogFile( Orbitals* _orbitals ) {

    return true;
}
     
/**
 * Converts the Gaussian data stored in the Orbitals object to GW input format
 */
bool GW::ConvertToGW( Orbitals* _orbitals ) {
    cerr << "Tried to convert to GW from package. ";
    throw std::runtime_error( "This does not make sense and should not have happened!");
}



}}
