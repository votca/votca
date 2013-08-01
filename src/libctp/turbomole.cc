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

#include "votca/ctp/turbomole.h"
#include "votca/ctp/segment.h"
#include <votca/tools/globals.h>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <stdio.h>
#include <iomanip>
#include <sys/stat.h>

using namespace std;

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
Turbomole::Turbomole( tools::Property *opt ) { 
               
    string key = "package";

    string _name = opt->get(key+".name").as<string> ();
    
    if ( _name != "turbomole" ) {
        cerr << "Tried to use " << _name << " package. ";
        throw std::runtime_error( "Package is not supported.");
    }
    
    _executable =       opt->get(key + ".executable").as<string> ();
    _options =          opt->get(key + ".options").as<string> ();
    _scratch_dir =      opt->get(key + ".scratch").as<string> ();
    _cleanup =          opt->get(key + ".cleanup").as<string> ();
    
    _write_guess = false;
    _get_charges = false;
    _get_self_energy = false;

    
};   
    

/**
 * Prepares the com file from a vector of segments
 * Appends a guess constructed from monomer orbitals if supplied
 */
bool Turbomole::WriteInputFile( vector<Segment* > segments, Orbitals* orbitals_guess )
{
    vector< Atom* > _atoms;
    vector< Atom* > ::iterator ait;
    vector< Segment* >::iterator sit;
    
    double nm2Bohr = 18.897259886;
    
    int qmatoms = 0;

    ofstream _com_file;
    
    string _com_file_name_full = _run_dir + "/" + _com_file_name;
    
    //cerr << "FILE NAME: " << _com_file_name_full << endl;
    
    _com_file.open ( _com_file_name_full.c_str() );

    // header 
    _com_file << "$coord" << endl;

    for (sit = segments.begin() ; sit != segments.end(); ++sit) {
        _atoms = (*sit)-> Atoms();

        for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

            if ((*ait)->HasQMPart() == false) { continue; }

            vec     pos = (*ait)->getQMPos();
            string  name = (*ait)->getElement();

            //fprintf(out, "%2s %4.7f %4.7f %4.7f \n"
            _com_file << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX()*nm2Bohr
                      << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY()*nm2Bohr
                      << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ()*nm2Bohr 
                      << setw(3) << name.c_str()   
                      << endl;
        }
    } 
    
    _com_file << "$end" << endl;
    _com_file.close();
    
    //cerr << _options << flush;
    
    LOG(logDEBUG,*_pLog) << "Preparing TURBOMOLE input " << flush;
    string _command;

    string _command_file_name = _run_dir + "/input";
    _com_file.open ( _command_file_name.c_str() );
    _com_file << "\n" << _options;
    _com_file.close();
        
    string _input_exe = "define";
    _command  = "cd " + _run_dir + "; " + _input_exe + " <  ./input" ;
    cerr << _command << flush;
    
    int i = system ( _command.c_str() );
    LOG(logDEBUG,*_pLog) << "Finished TURBOMOLE input" << flush;
    return true;    
    
}

bool Turbomole::WriteShellScript() {
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
 * Runs the Gaussian job. Returns 
 */
bool Turbomole::Run()
{

    LOG(logDEBUG,*_pLog) << "Running TURBOMOLE job" << flush;
    
    if (system(NULL)) {
        // if scratch is provided, run the shell script; 
        // otherwise run gaussian directly and rely on global variables 
        string _command;
        if ( _scratch_dir.size() != 0 ) {
            _command  = "cd " + _run_dir + "; tcsh " + _shell_file_name;
        }
        else {
            _command  = "cd " + _run_dir + "; " + _executable + " > " + _executable + ".log ";
        }
        
        int i = system ( _command.c_str() );
        LOG(logDEBUG,*_pLog) << "Finished TURBOMOLE job" << flush;
        return true;
    }
    else {
        LOG(logERROR,*_pLog) << _com_file_name << " failed to start" << flush; 
        return false;
    }
    



}

/**
 * Cleans up after the Gaussian job
 */
void Turbomole::CleanUp() {
    
    // cleaning up the generated files
    if ( _cleanup.size() != 0 ) {
        Tokenizer tok_cleanup(_cleanup, " \t\n");
        vector <string> _cleanup_info;
        tok_cleanup.ToVector(_cleanup_info);
        
        vector<string> ::iterator it;
               
        for (it = _cleanup_info.begin(); it != _cleanup_info.end(); ++it) {
            if ( *it == "com" ) {
                string file_name = _run_dir + "/" + _com_file_name;
                remove ( file_name.c_str() );
            }
            
            if ( *it == "sh" ) {
                string file_name = _run_dir + "/" + _shell_file_name;
                remove ( file_name.c_str() );
            }
            
            if ( *it == "log" ) {
                string file_name = _run_dir + "/" + _log_file_name;
                remove ( file_name.c_str() );
            }

           if ( *it == "chk" ) {
                string file_name = _run_dir + "/" + _chk_file_name;
                remove ( file_name.c_str() );
            }
            
            if ( *it == "fort.7" ) {
                string file_name = _run_dir + "/" + *it;
                remove ( file_name.c_str() );
            }            
        }
    }
    
}



/**
 * Reads in the MO coefficients from a GAUSSIAN fort.7 file
 */
bool Turbomole::ParseOrbitalsFile( Orbitals* _orbitals )
{
    std::map <int, std::vector<double> > _coefficients;
    std::map <int, double> _energies;
    
    std::string _line;
    unsigned _levels = 0;
    unsigned _level;
    unsigned _basis_size = 0;

    std::ifstream _input_file( _orb_file_name.c_str() );
    
    if (_input_file.fail()) {
        LOG( logERROR, *_pLog ) << "File " << _orb_file_name << " with molecular orbitals is not found " << flush;
        return false;
    } else {
        LOG(logDEBUG, *_pLog) << "Reading MOs from " << _orb_file_name << flush;
    }

    // skip lines with $ and #
    getline(_input_file, _line);
    getline(_input_file, _line);
    std::string::size_type hash_pos = _line.find("#");
    
    while ( hash_pos != std::string::npos ) {
           getline(_input_file, _line);
           hash_pos = _line.find("#");      
           //cout << _line;
    }
    
    //clog << endl << "Orbital file " << filename << " has " 
    //        << nrecords_in_line << " records per line, in D"
    //        << format << " format." << endl;

    while (_input_file) {

        // if a line has an equality sign, must be energy
        std::string::size_type energy_pos = _line.find("=");
        std::string::size_type dollar_pos = _line.find("$");

        if (energy_pos != std::string::npos  ) {

            vector<string> results;
            boost::trim( _line );
            
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                    boost::algorithm::token_compress_on); 
            //cout << results[1] << ":" << results[2] << ":" << results[3] << ":" << results[4] << endl;
            
            //exit(0);
            
            _level = boost::lexical_cast<int>(results.front());
            boost::replace_first(results[3], "D", "e");
            _energies[ _level ] = boost::lexical_cast<double>( results.back() );            
            _levels++;

        } else if ( dollar_pos == std::string::npos ) {
            
            while (_line.size() > 1) {
                string _coefficient;
                _coefficient.assign( _line, 0, 20 );
                boost::trim( _coefficient );
                boost::replace_first( _coefficient, "D", "e" );
                double coefficient = boost::lexical_cast<double>( _coefficient );
                _coefficients[ _level ].push_back( coefficient );
                _line.erase(0, 20);
            }
        }
        
        getline(_input_file, _line);
    }

    // some sanity checks
    LOG( logDEBUG, *_pLog ) << "Energy levels: " << _levels << flush;

    std::map< int, vector<double> >::iterator iter = _coefficients.begin();
    _basis_size = iter->second.size();

    for (iter = _coefficients.begin()++; iter != _coefficients.end(); iter++) {
        if (iter->second.size() != _basis_size) {
            LOG( logERROR, *_pLog ) << "Error reading " << _orb_file_name << ". Basis set size change from level to level." << flush;
            return false;
        }
    }
    
    LOG( logDEBUG, *_pLog ) << "Basis set size: " << _basis_size << flush;

    // copying information to the orbitals object
    _orbitals->_basis_set_size = _basis_size;
    _orbitals->_has_basis_set_size = true;
    _orbitals->_has_mo_coefficients = true;
    _orbitals->_has_mo_energies = true;
    
   // copying energies to a matrix  
   _orbitals->_mo_energies.resize( _levels );
   _level = 1;
   for(size_t i=0; i < _orbitals->_mo_energies.size(); i++) {
         _orbitals->_mo_energies[i] = _energies[ _level++ ];
   }
   
   // copying orbitals to the matrix
   (_orbitals->_mo_coefficients).resize( _levels, _basis_size );     
   for(size_t i = 0; i < _orbitals->_mo_coefficients.size1(); i++) {
      for(size_t j = 0 ; j < _orbitals->_mo_coefficients.size2(); j++) {
         _orbitals->_mo_coefficients(i,j) = _coefficients[i+1][j];
         //cout << i << " " << j << endl;
      }
   }

    
   //cout << _mo_energies << endl;   
   //cout << _mo_coefficients << endl; 
   
   // cleanup
   _coefficients.clear();
   _energies.clear();
   
     
   LOG(logDEBUG, *_pLog) << "Done reading MOs" << flush;

   return true;
}

bool Turbomole::CheckLogFile() {
    
    // check if the log file exists
    char ch;
    ifstream _input_file(_log_file_name.c_str());
    
    if (_input_file.fail()) {
        LOG(logERROR,*_pLog) << "Gaussian LOG is not found" << flush;
        return false;
    };

    _input_file.seekg(0,ios_base::end);   // go to the EOF
    
    // get empty lines and end of lines out of the way
    do {
        _input_file.seekg(-2,ios_base::cur);
        _input_file.get(ch);   
        //cout << "\nChar: " << ch << endl;
    } while ( ch == '\n' || ch == ' ' || (int)_input_file.tellg() == -1 );
 
    // get the beginning of the line or the file
    do {
       _input_file.seekg(-2,ios_base::cur);
       _input_file.get(ch);   
       //cout << "\nNext Char: " << ch << " TELL G " <<  (int)_input_file.tellg() << endl;
    } while ( ch != '\n' || (int)_input_file.tellg() == -1 );
            
    string _line;            
    getline(_input_file,_line);                      // Read the current line
    //cout << "\nResult: " << _line << '\n';     // Display it
    _input_file.close();
        
    std::string::size_type self_energy_pos = _line.find("Normal termination of Gaussian");
    if (self_energy_pos == std::string::npos) {
            LOG(logERROR,*_pLog) << "Gaussian LOG is incomplete" << flush;
            return false;      
    } else {
            //LOG(logDEBUG,*_pLog) << "Gaussian LOG is complete" << flush;
            return true;
    }
}

/**
 * Parses the Gaussian Log file and stores data in the Orbitals object 
 */
bool Turbomole::ParseLogFile( Orbitals* _orbitals ) {

    string _line;
    vector<string> results;
    bool _has_occupied_levels = false;
    bool _has_unoccupied_levels = false;
    bool _has_number_of_electrons = false;
    bool _has_basis_set_size = false;
    bool _has_overlap_matrix = false;
    bool _has_charges = false;
    bool _has_coordinates = false;
    bool _has_qm_energy = false;
    bool _has_self_energy = false;
    
    int _occupied_levels = 0;
    int _unoccupied_levels = 0;
    int _number_of_electrons = 0;
    int _basis_set_size = 0;

    
    LOG(logDEBUG,*_pLog) << "Parsing " << _log_file_name << flush;

    // check if LOG file is complete
    //if ( !CheckLogFile() ) return false;
    
    // Start parsing the file line by line
    ifstream _input_file(_log_file_name.c_str());
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
            string _oo = results.back();
            boost::trim( _oo );
            _number_of_electrons =  boost::lexical_cast<int>(_oo) ;
            _orbitals->_number_of_electrons = _number_of_electrons ;
            _orbitals->_has_number_of_electrons = true;
            LOG(logDEBUG,*_pLog) << "Alpha electrons: " << _number_of_electrons << flush ;
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
            string _bf = results.back();
            boost::trim( _bf );
            _basis_set_size = boost::lexical_cast<int>(_bf);
            _orbitals->_basis_set_size = _basis_set_size ;
            _orbitals->_has_basis_set_size = true;
            LOG(logDEBUG,*_pLog) << "Basis functions: " << _basis_set_size << flush;
        }

        /*
         * overlap matrix
         * stored after the *** Overlap *** line
         */
        std::string::size_type overlap_pos = _line.find("OVERLAP");
        if (overlap_pos != std::string::npos ) {

            // prepare the container
            _orbitals->_has_overlap = true;
            (_orbitals->_overlap).resize( _basis_set_size );
            
            _has_overlap_matrix = true;
            LOG(logDEBUG,*_pLog) << "Read the overlap matrix" << flush;
            
            // skip the next line with "----"
            getline(_input_file, _line);
            getline(_input_file, _line);
          
            overlap_pos = _line.find("--");
            
            int _i_index = 0;
            int _j_index = 0;
            
            while ( overlap_pos == std::string::npos ) {
                

                boost::trim( _line );
                
                vector<string> _row;
                boost::algorithm::split( _row , _line, boost::is_any_of(" "), boost::algorithm::token_compress_on);   
                int nfields =  _row.size();
                //cout << nfields << endl;
                
                for ( vector<string>::iterator it = _row.begin() ; it < _row.end() ; it++   ) {
                    
                    //cout << "  " << *it << endl;
                            
                    boost::trim( *it );
                    double _coefficient = boost::lexical_cast<double>( *it );
                    //cout << _i_index << ":" << _j_index << ":" << _coefficient << endl;
                    
                    _orbitals->_overlap( _i_index , _j_index ) = boost::lexical_cast<double>( _coefficient );
                    
                    _j_index++;
                    if ( _j_index > _i_index ) {
                        _j_index = 0;
                        _i_index++;
                    }
                }
                
                getline(_input_file, _line);
                overlap_pos = _line.find("--");

            } 
           

        } // end of the if "Overlap" found   

        
        /*
         *  Partial charges from the input file
         */
        std::string::size_type charge_pos = _line.find("Charges from ESP fit, RMS");
        
        if (charge_pos != std::string::npos && _get_charges ) {        
                LOG(logDEBUG,*_pLog) << "Getting charges" << flush;
                _has_charges = true;
                getline(_input_file, _line);
                getline(_input_file, _line);
                
                vector<string> _row;
                getline(_input_file, _line);
                boost::trim( _line );
                //cout << _line << endl;
                boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on); 
                int nfields =  _row.size();
                //cout << _row.size() << endl;
                
                while ( nfields == 3 ) {
                    int atom_id = boost::lexical_cast< int >( _row.at(0) );
                    int atom_number = boost::lexical_cast< int >( _row.at(0) );
                    string atom_type = _row.at(1);
                    double atom_charge = boost::lexical_cast< double >( _row.at(2) );
                    //if ( tools::globals::verbose ) cout << "... ... " << atom_id << " " << atom_type << " " << atom_charge << endl;
                    getline(_input_file, _line);
                    boost::trim( _line );
                    boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);  
                    nfields =  _row.size();
                    
                     if ( _orbitals->_has_atoms == false ) {
                         _orbitals->AddAtom( atom_type, 0, 0, 0, atom_charge );
                     } else {
                         QMAtom* pAtom = _orbitals->_atoms.at( atom_id - 1 );
                         pAtom->type = atom_type;
                         pAtom->charge = atom_charge;
                     }
                    
                }
                _orbitals->_has_atoms = true;
        }
        

         /*
         * Coordinates of the final configuration
         * stored in the archive at the end of the file
         */
         std::string::size_type coordinates_pos = _line.find("Test job not archived");
        
        if (coordinates_pos != std::string::npos) {
            LOG(logDEBUG,*_pLog) << "Getting the coordinates" << flush;
            _has_coordinates = true;
            string archive;
            while ( _line.size() != 0 ) {
                getline(_input_file, _line);
                boost::trim(_line);
                archive += _line;                
            }
            
            std::list<std::string> stringList;
            vector<string> results;
            boost::iter_split( stringList, archive, boost::first_finder("\\\\") );
             
            list<string>::iterator coord_block = stringList.begin();
            advance(coord_block, 3);
            
            vector<string> atom_block;
            boost::algorithm::split(atom_block, *coord_block, boost::is_any_of("\\"), boost::algorithm::token_compress_on);
            
            vector<string>::iterator atom_block_it;
            int aindex = 0;
            
            for(atom_block_it =  ++atom_block.begin(); atom_block_it != atom_block.end(); ++atom_block_it) {
                vector<string> atom;
                
                boost::algorithm::split(atom, *atom_block_it, boost::is_any_of(","), boost::algorithm::token_compress_on);
                string _atom_type = atom.front() ; 
                
                vector<string>::iterator it_atom;
                it_atom = atom.end();
                double _z =  boost::lexical_cast<double>( *(--it_atom) );
                double _y =  boost::lexical_cast<double>( *(--it_atom) );
                double _x =  boost::lexical_cast<double>( *(--it_atom) );
                
                if ( _orbitals->_has_atoms == false ) {
                        _orbitals->AddAtom( _atom_type, _x, _y, _z );
                } else {
                         QMAtom* pAtom = _orbitals->_atoms.at( aindex );
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
            boost::algorithm::split(energy, block[1], boost::is_any_of("="), boost::algorithm::token_compress_on);
            _orbitals->_qm_energy = _conv_Hrt_eV * boost::lexical_cast<double> ( energy[1] );
            
            LOG(logDEBUG, *_pLog) << "QM energy " << _orbitals->_qm_energy <<  flush;
                    
            _orbitals->_has_atoms = true;
            _orbitals->_has_qm_energy = true;

        }

         /*
         * Self-energy of external charges
         */
         std::string::size_type self_energy_pos = _line.find("Self energy of the charges");
        
        if (self_energy_pos != std::string::npos) {
            LOG(logDEBUG,*_pLog) << "Getting the self energy\n";  
            vector<string> block;
            vector<string> energy;
            boost::algorithm::split(block, _line, boost::is_any_of("="), boost::algorithm::token_compress_on);
            boost::algorithm::split(energy, block[1], boost::is_any_of("\t "), boost::algorithm::token_compress_on);
            
            _orbitals->_has_self_energy = true;
            _orbitals->_self_energy = _conv_Hrt_eV * boost::lexical_cast<double> ( energy[1] );
            
            LOG(logDEBUG, *_pLog) << "Self energy " << _orbitals->_self_energy <<  flush;

        }
        
        // check if all information has been accumulated and quit 
        if ( _has_number_of_electrons && 
             _has_basis_set_size && 
             _has_occupied_levels && 
             _has_unoccupied_levels &&
             _has_overlap_matrix &&
             _has_charges && 
             _has_self_energy
           ) break;
        
    } // end of reading the file line-by-line
   
    LOG(logDEBUG,*_pLog) << "Done parsing" << flush;
    return true;
}

string Turbomole::FortranFormat( const double &number ) {
    stringstream _ssnumber;
    if ( number >= 0) _ssnumber << " ";
    _ssnumber <<  setiosflags(ios::fixed) << setprecision(8) << std::scientific << number;
    std::string _snumber = _ssnumber.str(); 
    boost::replace_first(_snumber, "e", "D");
    return _snumber;
}
        



}}
