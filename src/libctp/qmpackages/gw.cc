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

     // ISOGWA file names
    string fileName = "system";

    _xyz_file_name = fileName + ".xyz";
    _input_file_name = "input"; // hard coded
    _log_file_name = fileName + ".gwbse"; 
    _orb_file_name = fileName + ".movecs" ;               

    string key = "package";
    string _name = options->get(key+".name").as<string> ();
    
    if ( _name != "gw" ) {
        cerr << "Tried to use " << _name << " package. ";
        throw std::runtime_error( "Wrong options file");
    }

    // set a default executable
    if ( options->exists(key + ".executable") ) {
        _executable =       options->get(key + ".executable").as<string> ();
    } else {
        _executable = "isogwa.x";
    } 
    
    // getting level ranges 
    _ranges = options->get(key + ".ranges").as<string> ();
    // now check validity, and get rpa, qp, and bse level ranges accordingly
    if ( _ranges == "default" ) {
        ;
    }else if ( _ranges == "factor" ) {
        // get factors
        _rpamaxfactor = options->get(key + ".rpamax").as<double> ();
        _qpminfactor  = options->get(key + ".qpmin").as<double> ();
        _qpmaxfactor  = options->get(key + ".qpmax").as<double> ();
        _bseminfactor = options->get(key + ".bsemin").as<double> ();
        _bsemaxfactor = options->get(key + ".bsemax").as<double> ();
    } else if ( _ranges == "explicit" ) {
        //get explicit numbers
        _rpamax = options->get(key + ".rpamax").as<unsigned int> ();
        _qpmin  = options->get(key + ".qpmin").as<unsigned int> ();
        _qpmax  = options->get(key + ".qpmax").as<unsigned int> ();
        _bsemin = options->get(key + ".bsemin").as<unsigned int> ();
        _bsemax = options->get(key + ".bsemax").as<unsigned int> ();
    } else {
        cerr << "\nSpecified range option " << _ranges << " invalid. ";
        throw std::runtime_error( "\nValid options are: default,factor,explicit");
    }
        
    _gwbasis =           options->get(key + ".gwbasis").as<string> ();

    
    //_charge =           options->get(key + ".charge").as<int> ();
    //_spin =             options->get(key + ".spin").as<int> ();
    _options =          options->get(key + ".options").as<string> ();
    _memory =           options->get(key + ".memory").as<string> ();
    _threads =          options->get(key + ".threads").as<int> ();
    _scratch_dir =      options->get(key + ".scratch").as<string> ();
    _cleanup =          options->get(key + ".cleanup").as<string> ();
    
}    

/**
 * Prepares the input file from a vector of segments
 * Appends a guess constructed from monomer orbitals if supplied
 */
bool GW::WriteInputFile( vector<Segment* > segments, Orbitals* orbitals_guess )
{
    vector< Atom* > _atoms;
    vector< Atom* > ::iterator ait;
    vector< Segment* >::iterator sit;
    
    int qmatoms = 0;

    // preparing file stream
    ofstream _gw_file;
    string _gw_file_name_full = _run_dir + "/" + _input_file_name;
    _gw_file.open ( _gw_file_name_full.c_str() );

    // load GW auxiliary basis set
    BasisSet gwbs;
    string basis_name( _gwbasis );
    gwbs.LoadBasisSet( basis_name );
    LOG(logDEBUG,*_pLog) << "Loaded Basis Set " << basis_name << flush;
    // header 
    /* if ( _options.size() ) _com_file <<  _options << endl ; */

    _gw_file << "start_dft votca " << endl;
    _gw_file << "lattice_structure sc " << endl;
    _gw_file << "lattice_constant_angstrom  1500.00 " << endl;
    _gw_file << "units angstrom " << endl;

    for (sit = segments.begin() ; sit != segments.end(); ++sit) {
        
        _atoms = (*sit)-> Atoms();

        for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

            if ((*ait)->HasQMPart() == false) { continue; }

            vec     pos = (*ait)->getQMPos();
            string  name = (*ait)->getElement();

            //fprintf(out, "%2s %4.7f %4.7f %4.7f \n"
            _gw_file << "atom " << setw(3) << name.c_str() 
                      << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX()*10
                      << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY()*10
                      << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ()*10 
                      << endl;
            
            // Now get the GW basis info for this element
            _gw_file << "basis_gwa ";
            Element* element = gwbs.getElement(name); 
            for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                        
                        Shell* shell = (*its);
                        // shell type, number primitives, scale factor
                        //_com_file << shell->getType() << " " << shell->getSize() << " " << shell->getScale() << endl;
                        // _el_file << shell->getType() << " " << shell->getSize() << " " << FortranFormat( shell->getScale() ) << endl;
                        for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                            GaussianPrimitive* gaussian = *itg;
                            //_com_file << gaussian->decay << " " << gaussian->contraction << endl;
                            _gw_file <<  gaussian->decay << " " << shell->getLmax() << " ";
                        }
                    }
            _gw_file << endl;
            
        }
    } 
    
    // some more defaults
    _gw_file << "dft_truncation_radius 500.0 " << endl;
    _gw_file << "screening   rpa_nosumrule_eps0  (0.d0, 1.d0)   1.001 " << endl;
    _gw_file << "tolerance_threecenter 1.d-9 " << endl;
    _gw_file << "scissors_operator_rigid_ryd  .3" << endl;

    // convert _rpamax if needed and write
    int _vbm = orbitals_guess->getNumberOfElectrons();
    if ( _ranges == "default" ){
        _rpamax = orbitals_guess->getNumberOfLevels();
    } else if ( _ranges == "factor" ) {
        _rpamax = _rpamaxfactor * _vbm;
    }
    _gw_file << "rpa_band_summation 1 " << _rpamax << endl;

    // VBM
    _gw_file << "mbpt_vbm " << _vbm << endl;
    
    // convert _qpmin and _qpmax if needed, and write
    if ( _ranges == "default" ){
        _qpmin = 1;
        _qpmax = 2 * _vbm;
    } else if ( _ranges == "factor" ) {
        _qpmin = _vbm - int( _qpminfactor * _vbm );
        _qpmax = _vbm + int( _qpmaxfactor * _vbm );
    }
    _gw_file << "gwa_bands_to_be_corrected " << _qpmin << " " << _qpmax << endl;

    // convert _bsemin and _bsemax if needed, and write
    if ( _ranges == "default" ){
        _bsemin = 1;
        _bsemax = 2 * _vbm;
    } else if ( _ranges == "factor" ) {
        _bsemin = _vbm - int( _bseminfactor * _vbm );
        _bsemax = _vbm + int( _bsemaxfactor * _vbm );
    }
    _gw_file << "bse_bands " << _bsemin << " " << _bsemax << endl;
    _gw_file << "bse_momentum_operator " << endl;
    if ( _options.size() ) _gw_file <<  _options << endl ; 
    _gw_file.close();
    
    return true;
}

/**
 * Runs the GW job. 
 */
bool GW::Run()
{

    LOG(logDEBUG,*_pLog) << "Running GWBSE job" << flush;
    // get the path to the shared folders with xml files
    char *gwbse_dir = getenv("GWBSEDIR");
    if( gwbse_dir == NULL) throw std::runtime_error("GWBSEDIR not set, cannot run isogwa.x.");
    
    if (system(NULL)) {
                
        // if threads is provided and > 1
        string _command;
        if ( _threads == 1 ) {
            _command  = "cd " + _run_dir + "; $GWBSEDIR/bin/" + _executable + " > " +  _log_file_name ;
        } else {
            _command  = "cd " + _run_dir + "; $GWBSEDIR/bin/" + _executable + " > " +  _log_file_name ;
            //_command  = "cd " + _run_dir + "; mpirun -np " +  boost::lexical_cast<string>(_threads) + " " + _executable + " " + _input_file_name + "> "+  _log_file_name ;
        }
        
        int i = system ( _command.c_str() );
        
        if ( CheckLogFile() ) {
            LOG(logDEBUG,*_pLog) << "Finished GWBSE job" << flush;
            return true;
        } else {
            LOG(logDEBUG,*_pLog) << "GWBSE failed" << flush;
        }
    }
    else {
        LOG(logERROR,*_pLog) << "GWBSE failed to start" << flush; 
        return false;
    }
    
    return true;

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
        LOG(logERROR,*_pLog) << "ISOGWA LOG is not found" << flush;
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
        
    std::string::size_type self_energy_pos = _line.find("INF total program runtime");
    if (self_energy_pos == std::string::npos) {
            LOG(logERROR,*_pLog) << "ISOGWA is incomplete" << flush;
            return false;      
    } else {
            //LOG(logDEBUG,*_pLog) << "Gaussian LOG is complete" << flush;
            return true;
    }
    
    _input_file.close();
    
}

/**
 * Parses the Gaussian Log file and stores data in the Orbitals object 
 */
bool GW::ParseLogFile( Orbitals* _orbitals ) {
    std::map <int, std::vector<double> > _coefficients;
    std::map <int, std::vector<double> > _energies;
    
    std::string _line;
    unsigned _levels = 0;
    unsigned _level;
    unsigned _basis_size = 0;
    unsigned _homoindex;
    unsigned _bse_lower;
    unsigned _bse_upper;
    unsigned _bse_states; 
    unsigned _bse_written_states; 
    vector<string> results; 
    
    // _store_qppert -> read from *.gwbse file
    string _log_file_name_full = _run_dir + "/" + _log_file_name;
    std::ifstream _input_file( _log_file_name_full.c_str() );
    
    if (_input_file.fail()) {
        LOG( logERROR, *_pLog ) << "File " << _log_file_name << " with GWBSE output not found " << flush;
        return false;
    } else {
        LOG(logDEBUG, *_pLog) << "Reading from " << _log_file_name << flush;
    }

    // go through file until 
//    getline(_input_file, _line);
  //  std::vector<string> strs;
   // boost::algorithm::split(strs, _line, boost::is_any_of("(D)"));
    //clog << strs.at(1) << endl;
   // int nrecords_in_line = boost::lexical_cast<int>(strs.at(1));
   // string format = strs.at(2);

    //clog << endl << "Orbital file " << filename << " has " 
    //        << nrecords_in_line << " records per line, in D"
    //        << format << " format." << endl;

    while (_input_file) {

        getline(_input_file, _line);
        // if a line has "GWA calculating correlation energies DONE!", relevant data follows
        std::string::size_type QPpert_pos = _line.find("GWA calculating correlation energies DONE!");

        if (QPpert_pos != std::string::npos) {
            // _orbitals->_has_QPpert = true;
            //_orbitals->_QP_levels_index.clear(); //make sure vector is empty
            // next line is header -> skip
            getline(_input_file, _line);
            _levels = 0;
            // now, data follows until "GWA setting up full QP Hamiltonian"
            bool _is_QPdata = true;
            while ( _is_QPdata ){
                getline(_input_file, _line);
//                std::string::size_type end_QPpert =
                if ( _line.find("GWA setting up full QP") != std::string::npos ){
                    _is_QPdata = false;
                } else {
                    int level;
                    boost::trim( _line );
                    _levels++;
                    boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                    boost::algorithm::token_compress_on); 
                    // get QP level index and save in orbitals
                    level = boost::lexical_cast<int>(results[2]);
                    //_orbitals->_QP_levels_index.push_back(boost::lexical_cast<int>(results[2]));
                    // get the five (DFT, S_x,S_c, V_xc, E_qp) energies in temporary map
                    for (size_t ite=0; ite<5; ite++) {
                        _energies[ _levels ].push_back(boost::lexical_cast<double>(results[ite+3])); // DFT energy
                    }
                }
            }
            // now copy energies in map into orbitals object matrix
            (_orbitals->_QPpert_energies).resize( _levels, 5 );
            for (size_t itl=0; itl < _levels; itl++){
               for (size_t ite=0; ite<5; ite++) {
                   _orbitals->_QPpert_energies(itl,ite) = _energies[itl+1][ite];
               }
            }
        }
        
        // find info about considered occupied and empty levels in BSE

        std::string::size_type homoindex_pos = _line.find("DFT HOMO index");
        //
        if (homoindex_pos != std::string::npos) {
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
            boost::algorithm::token_compress_on); 
            // get HOMO index
            _homoindex = boost::lexical_cast<int>(results.back());
            LOG( logDEBUG, *_pLog ) << "DFT HOMO index: " << _homoindex << flush;
        }
        
        std::string::size_type bse_lower_pos = _line.find("BSE lower");
        //
        if (bse_lower_pos != std::string::npos) {
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
            boost::algorithm::token_compress_on); 
            // get HOMO index
            _bse_lower = boost::lexical_cast<int>(results.back());
            LOG( logDEBUG, *_pLog ) << "BSE lowest level: " << _bse_lower << flush;
        }
        
        std::string::size_type bse_upper_pos = _line.find("BSE upper");
        //
        if (bse_upper_pos != std::string::npos) {
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
            boost::algorithm::token_compress_on); 
            // get HOMO index
            _bse_upper = boost::lexical_cast<int>(results.back());
            LOG( logDEBUG, *_pLog ) << "BSE upper level: " << _bse_upper << flush;
        }
        
        std::string::size_type bse_nstates_pos = _line.find("BSE number of states");
        //
        if (bse_nstates_pos != std::string::npos) {
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
            boost::algorithm::token_compress_on); 
            // get number os states
            _bse_states = boost::lexical_cast<int>(results.back());
            LOG( logDEBUG, *_pLog ) << "BSE number of states: " << _bse_states << flush;
        }
        
        std::string::size_type bse_wstates_pos = _line.find("BSE number of written states");
        //
        if (bse_wstates_pos != std::string::npos) {
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
            boost::algorithm::token_compress_on); 
            // get number of states written
            _bse_written_states = boost::lexical_cast<int>(results.back());
            LOG( logDEBUG, *_pLog ) << "BSE number of written states: " << _bse_written_states << flush;
        }
    }

    unsigned _bse_matrix_dim = (_homoindex - _bse_lower +1 ) * (_bse_upper - _homoindex );
    LOG( logDEBUG, *_pLog ) << "BSE matrix dimension: " << _bse_matrix_dim << flush;
    // build index array
    unsigned _i_bse = 0;
//    _orbitals->_BSE_levels_indices.resize( _bse_matrix_dim, 2 );
//    for ( unsigned _i_occupied = _bse_lower; _i_occupied <= _homoindex; _i_occupied++){
//        for (unsigned _i_empty = _homoindex+1; _i_empty <= _bse_upper; _i_empty++ ){
//            _orbitals->_BSE_levels_indices( _i_bse, 0 ) = _i_occupied;
//            _orbitals->_BSE_levels_indices( _i_bse, 1 ) = _i_empty;
//            _i_bse++;
 //        }
 //   }

    // some sanity checks
    LOG( logDEBUG, *_pLog ) << "QP perturbative levels: " << _levels << flush;
    _input_file.close();
    
    // now read from diag_QP.dat, if present
    string _diag_file_name = "diag_QP.dat";
    string _diag_file_name_full = _run_dir + "/" + _diag_file_name;
    std::ifstream _diag_file( _diag_file_name_full.c_str() );
    
    if (_diag_file.fail()) {
        LOG( logDEBUG, *_pLog ) << "File " << _diag_file_name << " not found, skipping " << flush;
    } else {
        //LOG(logDEBUG, *_pLog) << "Reading diagonalized QP data from " << _diag_file_name << flush;
        //vector<string> results;
        //_orbitals->_has_QPdiag = true;
        getline(_diag_file, _line); // first line is number of levels
        boost::trim( _line );
            
        boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
        boost::algorithm::token_compress_on); 
        _levels = boost::lexical_cast<int>(results.front());
        // sanity check
        //if ( _levels != _orbitals->_QP_levels_index.size() ){
        //   LOG( logERROR, *_pLog ) << "Unmatching number of QP levels in " << _log_file_name << "\n from *.gwbse: " << _orbitals->_QP_levels_index.size() << "; from diag_QP.dat: " << _levels << flush;
        //   return false;
        //}
        // now loop over all levels
        _orbitals->_QPdiag_coefficients.resize(_levels, _levels);
        for (size_t _i_level=0; _i_level < _levels; _i_level++ ){
            getline(_diag_file, _line); // first line has QP level index and energy
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
            boost::algorithm::token_compress_on); 
            _level = boost::lexical_cast<int>(results.front());
            _orbitals->_QPdiag_energies.push_back(boost::lexical_cast<double>(results.back()));
            // next _levels lines contain expansion coefficients
            for (size_t _i_coef = 0; _i_coef < _levels; _i_coef++ ){
                boost::trim( _line );
                boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                boost::algorithm::token_compress_on); 
                _orbitals->_QPdiag_coefficients(_i_level,_i_coef) = boost::lexical_cast<double>(results.back());
            }
        }
        _diag_file.close();
        LOG(logDEBUG, *_pLog) << "Diagonalized QPs: " << _levels << flush;
    }

    
    
    // now read from singlets.99, if present
    string _singlet_file_name = "singlets.99";
    string _singlet_file_name_full = _run_dir + "/" + _singlet_file_name;
    std::ifstream _singlet_file( _singlet_file_name_full.c_str() );
    
    if (_singlet_file.fail()) {
        LOG( logDEBUG, *_pLog ) << "File " << _singlet_file_name << " not found, skipping " << flush;
    } else {
        //LOG(logDEBUG, *_pLog) << "Reading singlet exciton data from " << _singlet_file_name << flush;
        //_orbitals->_has_BSE_singlets = true;
        _orbitals->_BSE_singlet_energies.resize( _bse_written_states );
        _orbitals->_BSE_singlet_coefficients.resize( _bse_written_states, _bse_matrix_dim );
        for ( size_t _i_state = 0; _i_state < _bse_written_states; _i_state++ ) {
            getline(_singlet_file, _line); // first line is energy
            boost::trim( _line );
            
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
            boost::algorithm::token_compress_on); 
            _orbitals->_BSE_singlet_energies[_i_state] =  boost::lexical_cast<double>(results.back());
            for (size_t _i_coef = 0; _i_coef < _bse_matrix_dim; _i_coef++ ){
                getline(_singlet_file, _line); // next lines are coefficients
                boost::trim( _line );
                boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                boost::algorithm::token_compress_on);
                // cout << _i_state << ":" << _i_coef << ":" << boost::lexical_cast<double>(results.back()) << endl;
                _orbitals->_BSE_singlet_coefficients(_i_state,_i_coef) = boost::lexical_cast<double>(results.back());
            }
        }
        LOG(logDEBUG, *_pLog) << "Singlet excitons: " << _bse_written_states << flush;
    }
    _singlet_file.close();

    
    // now read from triplets.99, if present
    string _triplet_file_name = "triplets.99";
    string _triplet_file_name_full = _run_dir + "/" + _triplet_file_name;
    std::ifstream _triplet_file( _triplet_file_name_full.c_str() );
    
    if (_triplet_file.fail()) {
        LOG( logDEBUG, *_pLog ) << "File " << _triplet_file_name << " not found, skipping " << flush;
    } else {
        //LOG(logDEBUG, *_pLog) << "Reading triplet exciton data from " << _triplet_file_name << flush;
        //_orbitals->_has_BSE_triplets = true;
        _orbitals->_BSE_triplet_energies.resize( _bse_written_states );
        _orbitals->_BSE_triplet_coefficients.resize( _bse_written_states, _bse_matrix_dim );
        for ( size_t _i_state = 0; _i_state < _bse_written_states; _i_state++ ) {
            getline(_triplet_file, _line); // first line is energy
            boost::trim( _line );
            
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
            boost::algorithm::token_compress_on); 
            _orbitals->_BSE_triplet_energies[_i_state] =  boost::lexical_cast<double>(results.back());
            for (size_t _i_coef = 0; _i_coef < _bse_matrix_dim; _i_coef++ ){
                getline(_triplet_file, _line); // next lines are coefficients
                boost::trim( _line );
                boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                boost::algorithm::token_compress_on);
                // cout << _i_state << ":" << _i_coef << ":" << boost::lexical_cast<double>(results.back()) << endl;
                _orbitals->_BSE_triplet_coefficients(_i_state,_i_coef) = boost::lexical_cast<double>(results.back());
            }
        }
        LOG(logDEBUG, *_pLog) << "Triplet excitons: " << _bse_written_states << flush;
    }
    _triplet_file.close();
        
   _coefficients.clear();
   _energies.clear();
  
     
   LOG(logDEBUG, *_pLog) << "Done reading GW-BSE data" << flush;

   return true;
    

    // _store_qpdiag -> read from diag_QP.dat
    // _store_singlets -> read from *.gwbse file (free interband, eh-interaction contrib.) and singlets.99 (coefficients)
    // _store_triplets -> read from triplets.99 (energies and coefficients), contributions?
    // return true;
}
     
/**
 * Converts the Gaussian data stored in the Orbitals object to GW input format
 */
bool GW::ConvertToGW( Orbitals* _orbitals ) {
    cerr << "Tried to convert to GW from package. ";
    throw std::runtime_error( "This does not make sense and should not have happened!");
}



}}
