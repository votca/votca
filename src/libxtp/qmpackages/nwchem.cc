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

#include "nwchem.h"
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

namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
void NWChem::Initialize( Property *options ) {

     // NWChem file names
    string fileName = "system";

    _xyz_file_name = fileName + ".xyz";
    _input_file_name = fileName + ".nw";
    _log_file_name = fileName + ".log"; 
    _shell_file_name = fileName + ".sh"; 
    _orb_file_name = fileName + ".movecs" ;               

    string key = "package";
    string _name = options->get(key+".name").as<string> ();
    
    if ( _name != "nwchem" ) {
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
    
     if (options->exists(key + ".outputVxc")) {
                _output_Vxc = options->get(key + ".outputVxc").as<bool> ();   
            }
             else _output_Vxc=false;
    // check whether options string contains vxc output, the _outputVxc is set to true
    std::string::size_type iop_pos = _options.find(" intermediate tXC matrix");
    if (iop_pos != std::string::npos) {
        if (_output_Vxc){
            cout << "=== You do not have to specify outputting Vxc twice. Next time remove the print ""intermediate tXC matrix"" part from your options string. Please continue"<< endl;
        }
        else{
            cout << "=== So you do not want to output Vxc but still put it in the options string? I will assume that you want to output Vxc, be more consistent next time. "<< endl;
        }
        _output_Vxc = true;   
    } 
    else if (_output_Vxc==true){
         _options=_options+"\n\ndft\nprint \"intermediate tXC matrix\"\nvectors input system.movecs\nnoscf\nend\ntask dft";
    }
    
    
    
    // check if the optimize keyword is present, if yes, read updated coords
    iop_pos = _options.find(" optimize");
    if (iop_pos != std::string::npos) {
        _is_optimization = true;
    } else
    {
        _is_optimization = false;
    }

    // check if the esp keyword is present, if yes, get the charges and save them
    iop_pos = _options.find(" esp");
    if (iop_pos != std::string::npos) {
        _get_charges = true;
    } else
    {
        _get_charges = false;
    }

    // check if the charge keyword is present, if yes, get the self energy and save it
    iop_pos = _options.find("set bq background");
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
 * Prepares the *.nw file from a vector of segments
 * Appends a guess constructed from monomer orbitals if supplied
 */
bool NWChem::WriteInputFile( vector<Segment* > segments, Orbitals* orbitals_guess )
{
    vector< Atom* > _atoms;
    vector< Atom* > ::iterator ait;
    vector< Segment* >::iterator sit;
    vector<string> results;
    //int qmatoms = 0;
    string temp_suffix = "/id";
    string scratch_dir_backup = _scratch_dir;
    ofstream _com_file;
    ofstream _crg_file;

    string _com_file_name_full = _run_dir + "/" + _input_file_name;
    string _crg_file_name_full = _run_dir + "/background.crg" ;
    
    _com_file.open ( _com_file_name_full.c_str() );
    // header 
   _com_file << "geometry noautoz noautosym" << endl;
 
   
   
   
   if ( !_write_charges ) {
   
    for (sit = segments.begin() ; sit != segments.end(); ++sit) {
        
        _atoms = (*sit)-> Atoms();
        temp_suffix = temp_suffix + "_" + boost::lexical_cast<string>( (*sit)->getId() );
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
    } 
    
    
    
   } else {
       
       // part of QM coordinates
       
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
       
                
                // part for the MM charge coordinates
                _crg_file.open ( _crg_file_name_full.c_str() );
                             for (it = qmatoms->begin(); it < qmatoms->end(); it++) {
                    if ((*it)->from_environment) {
                        boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
                        fmt % (*it)->x % (*it)->y % (*it)->z % (*it)->charge;
                        if ((*it)->charge != 0.0) _crg_file << fmt << endl;
                    }
                }

                _crg_file << endl;
                
                
                _crg_file.close();
       
   } // end different coordinate file data depending on QM/MM or pure QM

    _com_file << "end\n";
    
    // write charge of the molecule
    _com_file <<  "\ncharge " << _charge << "\n" ;
     
    // writing scratch_dir info
    if ( _scratch_dir != "" ) {
        
        LOG(logDEBUG,*_pLog) << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

        // boost::filesystem::create_directories( _scratch_dir + temp_suffix );
        string _temp( "scratch_dir " + _scratch_dir + temp_suffix + "\n" );
        _com_file << _temp ;
    }

    _com_file << _options << "\n";

    
    if ( _write_guess ) { 
        if ( orbitals_guess == NULL ) {
            throw std::runtime_error( "A guess for dimer orbitals has not been prepared.");
        } else {
            ofstream _orb_file;
            string _orb_file_name_full = _run_dir + "/" + _orb_file_name;
            // get name of temporary ascii file and open it
            boost::algorithm::split(results, _orb_file_name, boost::is_any_of("."), boost::algorithm::token_compress_on);
            string _orb_file_name_ascii = _run_dir + "/" + results.front() + ".mos";
            _orb_file.open ( _orb_file_name_ascii.c_str() );

            // header 
            _orb_file << "#generated by VOTCA\nbasisum\ngeomsum\n\nscf\nFri Sep 13 00:00:00 2013\nscf\n1\n\n8\nao basis\n1\n" << flush;   
            int _size_of_basis = (orbitals_guess->_mo_energies).size();
            // number of basis functions
            _orb_file << _size_of_basis << endl;
            // number of orbitals
            _orb_file << _size_of_basis << endl;

            
            vector<int> _sort_index;
            
            orbitals_guess->SortEnergies( &_sort_index );
            
            int level = 1;
            int ncolumns = 3;
            // write occupations as double in three columns
            // occupied levels
            int column = 1;
            for ( int i=0 ; i < orbitals_guess->_number_of_electrons ; i++ ){
                _orb_file <<  FortranFormat( 2.0 );
                if (column == ncolumns) { _orb_file << endl; column = 0; }
                column++;
            }
            // unoccupied levels
            for ( int i=orbitals_guess->_number_of_electrons ; i < _size_of_basis ; i++ ){
                _orb_file <<  FortranFormat( 0.0 );
                if (column == ncolumns) { _orb_file << endl; column = 0; }
                column++;
            }    
            // extra endl
            if ( column != 1 ) { _orb_file << endl; }
          
            // write all energies in same format
            column = 1;
            for ( vector< int > ::iterator soi = _sort_index.begin(); soi != _sort_index.end(); ++ soi ) {
                double _energy = (orbitals_guess->_mo_energies)[*soi] ;
                _orb_file <<  FortranFormat( _energy );
                if (column == ncolumns) { _orb_file << endl; column = 0; }
                column++;
            }
            if ( column != 1 ) _orb_file << endl;
            
            // write coefficients in same format
            for ( vector< int > ::iterator soi = _sort_index.begin(); soi != _sort_index.end(); ++ soi ) {
                ub::matrix_row< ub::matrix<double> > mr (orbitals_guess->_mo_coefficients, *soi);
                column = 1;
                for (unsigned j = 0; j < mr.size (); ++j) {
                        _orb_file <<  FortranFormat( mr[j] );
                        if (column == ncolumns) { _orb_file << endl; column = 0; }
                        column++;
                }
                
                level++;
                if ( column != 1 ) _orb_file << endl;
            } 
            _orb_file << " 0.0000   0.0000" << endl;
            _orb_file.close();
            
            // now convert this ascii file to binary
            string _command;
            _command = "cd "  + _run_dir + "; asc2mov 5000 system.mos system.movecs > convert.log";
            int i = std::system ( _command.c_str() );
            if ( i == 0 ) {
                LOG(logDEBUG, *_pLog) << "Converted MO file from ascii to binary"<< flush;
            } else  {
                LOG( logERROR, *_pLog ) << "Conversion of binary MO file to binary failed. " << flush;
                return false;
            }
            
            
        }
    
    
    }   

    _com_file << endl;
    _com_file.close();
    
    
    
      // and now generate a shell script to run both jobs, if neccessary
      LOG(logDEBUG, *_pLog) << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;

      _scratch_dir = scratch_dir_backup + temp_suffix;
            
            //boost::filesystem::create_directories(_scratch_dir + temp_suffix);
            //string _temp("scratch_dir " + _scratch_dir + temp_suffix + "\n");
            //_com_file << _temp;
            WriteShellScript();
            _scratch_dir = scratch_dir_backup;
    
    
    
    
    
    
    return true;
    
}

bool NWChem::WriteShellScript() {
    ofstream _shell_file;
    
    string _shell_file_name_full = _run_dir + "/" + _shell_file_name;
            
    _shell_file.open ( _shell_file_name_full.c_str() );

    _shell_file << "#!/bin/bash" << endl ;
    _shell_file << "mkdir -p " << _scratch_dir << endl;
    
    if ( _threads == 1 ){ 
        _shell_file << _executable << " " << _input_file_name << " > " <<  _log_file_name  << " 2> run.error" << endl;    
    } else {
        _shell_file << "mpirun -np " << boost::lexical_cast<string>(_threads) << " " << _executable << " " << _input_file_name << " > " <<  _log_file_name << " 2> run.error" << endl;    
    }
    _shell_file.close();
    
    return true;   
}

/**
 * Runs the NWChem job. 
 */
bool NWChem::Run()
{

    LOG(logDEBUG,*_pLog) << "Running NWChem job" << flush;
    
    if (std::system(NULL)) {
        
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
            //_command  = "cd " + _run_dir + "; mkdir -p " + _scratch_dir + "; " + _executable + " " + _input_file_name + "> " +  _log_file_name ;
                        _command = "cd " + _run_dir + "; sh " + _shell_file_name;
        } else {
            //_command  = "cd " + _run_dir + "; mkdir -p " + _scratch_dir + ";  mpirun -np " +  boost::lexical_cast<string>(_threads) + " " + _executable + " " + _input_file_name + "> "+  _log_file_name ;
            _command = "cd " + _run_dir + "; sh " + _shell_file_name;
        }
        
        //int i = std::system ( _command.c_str() );
        std::system ( _command.c_str() );
        
        if ( CheckLogFile() ) {
            LOG(logDEBUG,*_pLog) << "Finished NWChem job" << flush;
            return true;
        } else {
            LOG(logDEBUG,*_pLog) << "NWChem job failed" << flush;
            return false;
        }
    }
    else {
        LOG(logERROR,*_pLog) << _input_file_name << " failed to start" << flush; 
        return false;
    }
    
    return true;
}

/**
 * Cleans up after the NWChem job
 */
void NWChem::CleanUp() {
    
    // cleaning up the generated files
    if ( _cleanup.size() != 0 ) {
        Tokenizer tok_cleanup(_cleanup, ",");
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
                string file_name = _run_dir + "/" + _orb_file_name;
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
bool NWChem::ParseOrbitalsFile( Orbitals* _orbitals )
{
    std::map <int, std::vector<double> > _coefficients;
    std::map <int, double> _energies;
    std::map <int, double> _occ;
    
    std::string _line;
    unsigned _levels = 0;
    //unsigned _level;
    unsigned _basis_size = 0;
    int _number_of_electrons = 0;
    
    /* maybe we DO need to convert from fortran binary to ASCII first to avoid
     compiler-dependent issues */
    string _orb_file_name_bin = _run_dir + "/" + _orb_file_name;
    string _command;
    _command = "cd "  + _run_dir + "; mov2asc 10000 system.movecs system.mos > convert.log";
    int i = std::system ( _command.c_str() );
    if ( i == 0 ) {
        LOG(logDEBUG, *_pLog) << "Converted MO file from binary to ascii"<< flush;
    } else  {
        LOG( logERROR, *_pLog ) << "Conversion of binary MO file to ascii failed. " << flush;
        return false;
    }
        
    // opening the ascii MO file   
    string _orb_file_name_full = _run_dir + "/" + "system.mos" ;
    std::ifstream _input_file( _orb_file_name_full.c_str() );
    
    if (_input_file.fail()) {
        LOG( logERROR, *_pLog ) << "File " << _orb_file_name_full << " with molecular orbitals is not found " << flush;
        return false;
    } else {
        LOG(logDEBUG, *_pLog) << "Reading MOs from " << _orb_file_name_full << flush;
    }
    
    // the first 12 lines are garbage info
    for ( i=1; i<13; i++) {
       getline(_input_file, _line);
    }
    // next line has basis set size
    _input_file >> _basis_size ;
    LOG( logDEBUG, *_pLog ) << "Basis set size: " << _basis_size << flush;


    // next line has number of stored MOs
    _input_file >> _levels ;
    LOG( logDEBUG, *_pLog ) << "Energy levels: " << _levels << flush;

    /* next lines contain information about occupation of the MOs
     *  - each line has 3 numbers 
     *  - from _occ we can determine the number of electrons/2 */
    int _n_lines = (( _levels - 1 ) / 3);
    int _n_rest  = _levels - 3*_n_lines;
    // read in the data
    int _imo = 0;
    for ( i=1; i<=_n_lines; i++ ) {
        for ( int j=0; j<3; j++ ){
            _input_file >> _occ[ _imo ];
            if ( _occ[ _imo ] == 2.0 ) {
                _number_of_electrons++;
            }
            _imo++;
        }
    }
    if ( _n_rest != 0 ) {
        for ( i=0; i<_n_rest; i++ ){
             _input_file >> _occ[ _imo ];
            _imo++;
        }
    }
    LOG(logDEBUG,*_pLog) << "Alpha electrons: " << _number_of_electrons << flush ;

    int _occupied_levels = _number_of_electrons;
    int _unoccupied_levels = _levels - _occupied_levels;
    LOG(logDEBUG,*_pLog) << "Occupied levels: " << _occupied_levels << flush;
    LOG(logDEBUG,*_pLog) << "Unoccupied levels: " << _unoccupied_levels << flush;  
    
    // reset index and read MO energies the same way
    _imo = 0;
    for ( i=1; i<=_n_lines; i++ ) {
        for ( int j=0; j<3; j++ ){
            _input_file >> _energies[ _imo ];
            _imo++;
        }
    }
    if ( _n_rest != 0 ) {
        for ( i=0; i<_n_rest; i++ ){
             _input_file >> _energies[ _imo ];
            _imo++;
        }
    }


    // Now, the same for the coefficients
    double coef;
    for ( unsigned _imo=0; _imo < _levels ; _imo++ ){
        for ( i=1; i<=_n_lines; i++ ) {
            for ( int j=0; j<3; j++ ){
                _input_file >> coef;
                _coefficients[ _imo ].push_back( coef );
            }
        }
        if ( _n_rest != 0 ) {
           for ( i=0; i<_n_rest; i++ ){
                _input_file >> coef;
                _coefficients[ _imo ].push_back( coef );
           }
        }
    }
    
                   

  
    // copying information to the orbitals object
    _orbitals->setBasisSetSize(  _basis_size );
    //_orbitals->_has_basis_set_size = true;
    _orbitals->setNumberOfElectrons( _number_of_electrons );
    // _orbitals->_has_number_of_electrons = true;
    // _orbitals->_has_mo_coefficients = true;
    // _orbitals->_has_mo_energies = true;
    //_orbitals->_has_occupied_levels = true;
    //_orbitals->_has_unoccupied_levels = true;
    //_orbitals->_occupied_levels = _occupied_levels;
    //_orbitals->_unoccupied_levels = _unoccupied_levels;
    _orbitals->setNumberOfLevels( _occupied_levels , _unoccupied_levels );
    
   // copying energies to a matrix  
   _orbitals->_mo_energies.resize( _levels );
   //_level = 1;
   for(size_t i=0; i < _orbitals->_mo_energies.size(); i++) {
         _orbitals->_mo_energies[i] = _energies[ i ];
   }

   
   // copying orbitals to the matrix
   (_orbitals->_mo_coefficients).resize( _levels, _basis_size );     
   for(size_t i = 0; i < _orbitals->_mo_coefficients.size1(); i++) {
      for(size_t j = 0 ; j < _orbitals->_mo_coefficients.size2(); j++) {
         _orbitals->_mo_coefficients(i,j) = _coefficients[i][j];
         //cout << i << " " << j << endl;
      }
   }

   
  
    
   //cout << _mo_energies << endl;   
   //cout << _mo_coefficients << endl; 
   
   // cleanup
   _coefficients.clear();
   _energies.clear();
   _occ.clear();
   
   // when all is done, we can trash the ascii file
   string file_name = _run_dir + "/system.mos";
   remove ( file_name.c_str() );
   
   
   
   LOG(logDEBUG, *_pLog) << "Done reading MOs" << flush;

   return true;
}

bool NWChem::CheckLogFile() {
    
    // check if the log file exists
    char ch;
    ifstream _input_file((_run_dir + "/" + _log_file_name).c_str());
    
    if (_input_file.fail()) {
        LOG(logERROR,*_pLog) << "NWChem LOG is not found" << flush;
        return false;
    };

   if (_input_file.peek() == std::ifstream::traits_type::eof()) {
       LOG(logERROR,*_pLog) << "NWChem run failed. Check OpenMPI version!" << flush;
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
       //cout << "\nResult: " << _line << '\n';     // Display it
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
    return true;    
}

/**
 * Parses the Gaussian Log file and stores data in the Orbitals object 
 */
bool NWChem::ParseLogFile( Orbitals* _orbitals ) {

    static const double _conv_Hrt_eV = 27.21138386;

    string _line;
    vector<string> results;

    bool _has_overlap_matrix = false;
    bool _has_charges = false;
    //bool _has_coordinates = false;
    //bool _has_vxc_matrix = false;
    bool _has_qm_energy = false;
    bool _has_self_energy = false;
    bool _has_basis_set_size = false;
    
    bool _found_optimization = false;
    int _basis_set_size = 0;
     
    LOG(logDEBUG,*_pLog) << "Parsing " << _log_file_name << flush;
       // return true;
    string _log_file_name_full =  _run_dir + "/" + _log_file_name;
    // check if LOG file is complete
    if ( !CheckLogFile() ) return false;
    
    // save qmpackage name
    // _orbitals->_has_qm_package = true;
    _orbitals->setQMpackage("nwchem");
    
    // set _found_optimization to true if this is a run without optimization
    if ( !_is_optimization ) {
        _found_optimization = true;
    }
    
    // Start parsing the file line by line
    ifstream _input_file(_log_file_name_full.c_str());
    while (_input_file) {

        getline(_input_file, _line);
        boost::trim(_line);

        /*
         * basis set size (is only required for overlap matrix reading, rest is
         * in orbitals file and could be skipped
         */
        std::string::size_type basis_pos = _line.find("number of functions");
        if (basis_pos != std::string::npos) {
            //cout << _line << endl;
            boost::algorithm::split(results, _line, boost::is_any_of(":"), boost::algorithm::token_compress_on);
            _has_basis_set_size = true;
            string _bf = results.back();
            boost::trim( _bf );
            _basis_set_size = boost::lexical_cast<int>(_bf);
            _orbitals->setBasisSetSize( _basis_set_size );
            LOG(logDEBUG,*_pLog) << "Basis functions: " << _basis_set_size << flush;
        }

      /*
       * Total DFT energy
       */
       std::string::size_type energy_pos = _line.find("Total DFT energy");
       if (energy_pos != std::string::npos) {
            //cout << _line << endl;
            boost::algorithm::split(results, _line, boost::is_any_of("="), boost::algorithm::token_compress_on);
            string _energy = results.back();
            boost::trim( _energy );
            _orbitals->setQMEnergy ( _conv_Hrt_eV * boost::lexical_cast<double>(_energy) );
            LOG(logDEBUG, *_pLog) << "QM energy " << _orbitals->getQMEnergy() <<  flush;
            _has_qm_energy = true;
            // _orbitals->_has_qm_energy = true;
            
        }
 
       
       /*
        * Vxc matrix
        * stored after the global array: g vxc
        */
       
       std::string::size_type vxc_pos = _line.find("global array: g vxc");
       if (vxc_pos != std::string::npos ) {
       
       
                // prepare the container
                // _orbitals->_has_vxc = true;
                ub::symmetric_matrix<double>& _vxc = _orbitals->AOVxc();
                _vxc.resize(_basis_set_size);


                //_has_vxc_matrix = true;
                          vector<int> _j_indeces;
            
           int _n_blocks = 1 + (( _basis_set_size - 1 ) / 6);
           //cout << _n_blocks;
            
           for (int _block = 0; _block < _n_blocks; _block++ ) {
                // first line is garbage
                getline(_input_file, _line); 
                // second line gives the j index in the matrix
                getline(_input_file, _line); 
                boost::tokenizer<> tok( _line );
                
                /// COMPILATION IS BROKEN DUE TO A BUG IN BOOST 1.53
                std::transform( tok.begin(), tok.end(), std::back_inserter( _j_indeces ), &boost::lexical_cast<int,std::string> );

                // third line is garbage again
                getline(_input_file, _line);

                // read the block of max _basis_size lines + the following header
                for (int i = 0; i < _basis_set_size; i++ ) {
                    getline (_input_file, _line); 
                    
                    // split the line on the i index and the rest
                    vector<string> _row;
                    boost::trim( _line );
                    boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on); 
                   
                            
                    int _i_index = boost::lexical_cast<int>(  _row.front()  ); 
                    _row.erase( _row.begin() );
                    
                    std::vector<int>::iterator _j_iter = _j_indeces.begin();
                    
                    for (std::vector<string>::iterator iter = _row.begin()++; iter != _row.end(); iter++) {
                        string  _coefficient = *iter;
                       
                        int _j_index = *_j_iter;                                
                        // _orbitals->_overlap( _i_index-1 , _j_index-1 ) = boost::lexical_cast<double>( _coefficient );
                        _vxc(_i_index-1 , _j_index-1) = boost::lexical_cast<double>( _coefficient );
                        _j_iter++;
                        
                    }
 
                    
                }
                
                // clear the index for the next block
                _j_indeces.clear();        
            } // end of the blocks
            LOG(logDEBUG,*_pLog) << "Read the Vxc matrix" << flush;
                
       }
       
                       /* Check for ScaHFX = factor of HF exchange included in functional */
                std::string::size_type HFX_pos = _line.find("Hartree-Fock (Exact) Exchange");
                if (HFX_pos != std::string::npos) {
                    boost::algorithm::split(results, _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);
                    double _ScaHFX = boost::lexical_cast<double>(results.back());
                    _orbitals->setScaHFX(_ScaHFX);
                    LOG(logDEBUG, *_pLog) << "DFT with " << _ScaHFX << " of HF exchange!" << flush;
                }
       
       
       
       
        /*
         * overlap matrix
         * stored after the global array: Temp Over line
         */
       std::string::size_type overlap_pos = _line.find("global array: Temp Over");
       if (overlap_pos != std::string::npos ) {
            
           // prepare the container
           // _orbitals->_has_overlap = true;
           (_orbitals->_overlap).resize( _basis_set_size );
            
           _has_overlap_matrix = true;
           vector<int> _j_indeces;
            
           int _n_blocks = 1 + (( _basis_set_size - 1 ) / 6);
           //cout << _n_blocks;
            
           for (int _block = 0; _block < _n_blocks; _block++ ) {
                // first line is garbage
                getline(_input_file, _line); 
                // second line gives the j index in the matrix
                getline(_input_file, _line); 
                boost::tokenizer<> tok( _line );
                
                /// COMPILATION IS BROKEN DUE TO A BUG IN BOOST 1.53
                std::transform( tok.begin(), tok.end(), std::back_inserter( _j_indeces ), &boost::lexical_cast<int,std::string> );

                // third line is garbage again
                getline(_input_file, _line);

                // read the block of max _basis_size lines + the following header
                for (int i = 0; i < _basis_set_size; i++ ) {
                    getline (_input_file, _line); 
                    
                    // split the line on the i index and the rest
                    vector<string> _row;
                    boost::trim( _line );
                    boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on); 
                   
                            
                    int _i_index = boost::lexical_cast<int>(  _row.front()  ); 
                    _row.erase( _row.begin() );
                    
                    std::vector<int>::iterator _j_iter = _j_indeces.begin();
                    
                    for (std::vector<string>::iterator iter = _row.begin()++; iter != _row.end(); iter++) {
                        string  _coefficient = *iter;
                       
                        int _j_index = *_j_iter;                                
                        _orbitals->_overlap( _i_index-1 , _j_index-1 ) = boost::lexical_cast<double>( _coefficient );
                        _j_iter++;
                        
                    }
 
                    
                }
                
                // clear the index for the next block
                _j_indeces.clear();        
            } // end of the blocks
            LOG(logDEBUG,*_pLog) << "Read the overlap matrix" << flush;
        } // end of the if "Overlap" found   

        
        /*
         *  Partial charges from the input file
         */
        std::string::size_type charge_pos = _line.find("ESP");
        

        if (charge_pos != std::string::npos && _get_charges ) {        
                LOG(logDEBUG,*_pLog) << "Getting charges" << flush;
                _has_charges = true;
                // two empty lines
                getline(_input_file, _line);
                getline(_input_file, _line);
                
                // now starts the data in format
                // _id type x y z q
                vector<string> _row;
                getline(_input_file, _line);
                boost::trim( _line );
                //cout << _line << endl;
                boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on); 
                int nfields =  _row.size();
                //cout << _row.size() << endl;
                
                while ( nfields == 6 ) {
                    int atom_id = boost::lexical_cast< int >( _row.at(0) );
                    //int atom_number = boost::lexical_cast< int >( _row.at(0) );
                    string atom_type = _row.at(1);
                    double atom_charge = boost::lexical_cast< double >( _row.at(5) );
                    //if ( tools::globals::verbose ) cout << "... ... " << atom_id << " " << atom_type << " " << atom_charge << endl;
                    getline(_input_file, _line);
                    boost::trim( _line );
                    boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);  
                    nfields =  _row.size();
                    
                     if ( _orbitals->hasQMAtoms() == false ) {
                         _orbitals->AddAtom( atom_type, 0, 0, 0, atom_charge );
                     } else {
                         QMAtom* pAtom = _orbitals->_atoms.at( atom_id - 1 );
                         pAtom->type = atom_type;
                         pAtom->charge = atom_charge;
                     }
                    
                }
                //_orbitals->_has_atoms = true;
        }
        

         /*
         * Coordinates of the final configuration
         * depending on whether it is an optimization or not
         */
        
        
        if ( _is_optimization ){
                std::string::size_type optimize_pos = _line.find("Optimization converged");
                if (optimize_pos != std::string::npos) {
                        _found_optimization = true;
                }
        }
        
        std::string::size_type coordinates_pos = _line.find("Output coordinates");
        
        if ( _found_optimization && coordinates_pos != std::string::npos) {
            LOG(logDEBUG,*_pLog) << "Getting the coordinates" << flush;
            
            //_has_coordinates = true;
            bool _has_QMAtoms = _orbitals->hasQMAtoms();

            // three garbage lines
            getline(_input_file, _line);
            getline(_input_file, _line);
            getline(_input_file, _line);
            // now starts the data in format
            // _id type Qnuc x y z 
            vector<string> _row;
            getline(_input_file, _line);
            boost::trim( _line );

            boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on); 
            int nfields =  _row.size();

                
            while ( nfields == 6 ) {
                int atom_id = boost::lexical_cast< int >( _row.at(0) );
                //int atom_number = boost::lexical_cast< int >( _row.at(0) );
                string _atom_type = _row.at(1);
                double _x =  boost::lexical_cast<double>( _row.at(3) );
                double _y =  boost::lexical_cast<double>( _row.at(4) );
                double _z =  boost::lexical_cast<double>( _row.at(5) );
                //if ( tools::globals::verbose ) cout << "... ... " << atom_id << " " << atom_type << " " << atom_charge << endl;
                getline(_input_file, _line);
                boost::trim( _line );
                boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on);  
                nfields =  _row.size();

                if ( _has_QMAtoms == false ) {
                    _orbitals->AddAtom( _atom_type, _x, _y, _z );
                } else {
                                       
                    QMAtom* pAtom = _orbitals->_atoms.at( atom_id - 1 );
                    pAtom->type = _atom_type;
                    pAtom->x = _x;
                    pAtom->y = _y;
                    pAtom->z = _z;
                }
                    
            }
      
            // _orbitals->_has_atoms = true;

        }
        
         /*
         * TODO Self-energy of external charges
         */
         std::string::size_type self_energy_pos = _line.find("Self energy of the charges");
        
        if (self_energy_pos != std::string::npos) {
            LOG(logDEBUG,*_pLog) << "Getting the self energy\n";  
            vector<string> block;
            vector<string> energy;
            boost::algorithm::split(block, _line, boost::is_any_of("="), boost::algorithm::token_compress_on);
            boost::algorithm::split(energy, block[1], boost::is_any_of("\t "), boost::algorithm::token_compress_on);
            
            // _orbitals->_has_self_energy = true;
            _orbitals->setSelfEnergy( _conv_Hrt_eV * boost::lexical_cast<double> ( energy[1] ) );
            
            LOG(logDEBUG, *_pLog) << "Self energy " << _orbitals->getSelfEnergy() <<  flush;

        }
        
        // check if all information has been accumulated and quit 
        if ( _has_basis_set_size && 
             _has_overlap_matrix &&
             _has_charges && 
             _has_qm_energy &&
             _has_self_energy
           ) break;
        
    } // end of reading the file line-by-line
   
    LOG(logDEBUG,*_pLog) << "Done parsing" << flush;
    return true;
}



/**
 * Converts the NWChem data stored in the Orbitals object to GW input format
 */
bool NWChem::ConvertToGW( Orbitals* _orbitals ) {
    cerr << "Tried to convert to GW from NWChem package. ";
    throw std::runtime_error( "Conversion not implemented yet!");
}


string NWChem::FortranFormat( const double &number ) {
    stringstream _ssnumber;
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
