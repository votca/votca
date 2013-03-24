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

#include "votca/ctp/orbitals.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <iterator>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;

namespace votca { namespace ctp {

Orbitals::Orbitals() { 
    
    _basis_set_size = 0;
    _occupied_levels = 0;
    _unoccupied_levels = 0;
    _electrons = 0;
    
    _has_basis_set_size = false;
    _has_occupied_levels = false;
    _has_unoccupied_levels = false;
    _has_electrons = false;   
    _has_degeneracy = false;
    _has_mo_coefficients = false;
    _has_mo_energies = false;
    _has_overlap = false;
    
};   
    
Orbitals::~Orbitals() { 
    _mo_energies.clear();
    _mo_coefficients.clear();
    _overlap.clear();
};   

/**
 * Reads in the MO coefficients from a GAUSSIAN fcheck file
 */
bool Orbitals::ReadOrbitalsGaussian( const char * filename )
{
    std::map <int, std::vector<double> > _coefficients;
    std::map <int, double> _energies;
    
    std::string _line;
    unsigned _levels = 0;
    unsigned _level;
    unsigned _basis_size = 0;

    std::ifstream _input_file(filename);
    if (_input_file.fail()) {
        cerr << endl << "File " << filename << " with molecular orbitals is not found " << endl;
        return 1;
    } else {
        cout << endl << "... reading molecular orbitals from " << filename << endl;
    }

    // number of coefficients per line is  in the first line of the file (5D15.8)
    getline(_input_file, _line);
    std::vector<string> strs;
    boost::algorithm::split(strs, _line, boost::is_any_of("(D)"));
    //cout << strs.at(1) << endl;
    int nrecords_in_line = boost::lexical_cast<int>(strs.at(1));
    string format = strs.at(2);

    //cout << endl << "Orbital file " << filename << " has " 
    //        << nrecords_in_line << " records per line, in D"
    //        << format << " format." << endl;

    while (_input_file) {

        getline(_input_file, _line);
        // if a line has an equality sign, must be energy
        std::string::size_type energy_pos = _line.find("=");

        if (energy_pos != std::string::npos) {

            vector<string> results;
            boost::trim( _line );
            
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                    boost::algorithm::token_compress_on); 
            //cout << results[1] << ":" << results[2] << ":" << results[3] << ":" << results[4] << endl;
            
            _level = boost::lexical_cast<int>(results.front());
            boost::replace_first(results.back(), "D", "e");
            _energies[ _level ] = boost::lexical_cast<double>( results.back() );            
            _levels++;

        } else {
            
            while (_line.size() > 1) {
                string _coefficient;
                _coefficient.assign( _line, 0, 15 );
                boost::trim( _coefficient );
                boost::replace_first( _coefficient, "D", "e" );
                double coefficient = boost::lexical_cast<double>( _coefficient );
                _coefficients[ _level ].push_back( coefficient );
                _line.erase(0, 15);
            }
        }
    }

    // some sanity checks
    cout << "... ... energy levels: " << _levels << endl;

    std::map< int, vector<double> >::iterator iter = _coefficients.begin();
    _basis_size = iter->second.size();

    for (iter = _coefficients.begin()++; iter != _coefficients.end(); iter++) {
        if (iter->second.size() != _basis_size) {
            cerr << "Error reading " << filename << ". Basis set size change from level to level.";

        }
    }
    cout << "... ... basis set size: " << _basis_size << endl;

    _basis_set_size = _basis_size;
    _has_basis_set_size = true;
    
    
   // copying energies to a matrix
   _mo_energies.resize( _levels );
   _level = 1;
   for(size_t i=0; i<_mo_energies.size(); i++) {
         _mo_energies[i] = _energies[ _level++ ];
   }
   
   // copying orbitals to the matrix
   _mo_coefficients.resize( _levels, _basis_size );     
   for(size_t i = 0; i < _mo_coefficients.size1(); i++) {
      for(size_t j = 0 ; j<_mo_coefficients.size2(); j++) {
         _mo_coefficients(i,j) = _coefficients[i+1][j];
         //cout << i << " " << j << endl;
      }
   }
   
   //cout << _mo_energies << endl;   
   //cout << _mo_coefficients << endl; 
   
   // cleanup
   _coefficients.clear();
   _energies.clear();
   
     
   cout << "... done reading orbital files from " << filename << endl;
   _has_mo_coefficients = true;
   _has_mo_energies = true;
   return 0;
}


/**
 * Reads in the Orbital Overlap matrix from a GAUSSIAN log file
 */
bool Orbitals::ReadOverlapGaussian( const char * filename )
{
    string _line;
    unsigned _basis_size = 0;
      
    ifstream _input_file(filename);
    if (_input_file.fail()) {
        cerr << endl << "File " << filename << " with overlap matrix is not found " << endl;
        return 1;
    } else {
        cout << endl << "... reading the overlap matrix from " << filename << endl;
    }
    
    while (_input_file ) {

        getline(_input_file, _line);
        // if a line has an equality sign, must be energy
        std::string::size_type overlap_pos = _line.find("Overlap");
        std::string::size_type nbasis_pos = _line.find("NBasis");
 
        if (nbasis_pos != std::string::npos && _basis_size == 0 ) {
          
            vector<string> results;
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                    boost::algorithm::token_compress_on); 
            //cout << results[1] << ":" << results[2] << ":" << results[3] << ":" << results[4] << endl;
            
            _basis_size = boost::lexical_cast<int>(results[1]);
            cout << "... ... number of basis functions: " << _basis_size << endl;
            // preparing the matrix container
            _overlap.resize( _basis_size );  
 
        }
        
        // get out of this loop if we already have the overlap matrix
        if ( _has_overlap ) { break; }
                    
        if (overlap_pos != std::string::npos ) {
            
            _has_overlap = true;
            //cout << "Found the overlap matrix!" << endl;   
            vector<int> _j_indeces;
            
            int _n_blocks = 1 + (( _basis_size - 1 ) / 5);
            //cout << _n_blocks;
            
            getline(_input_file, _line); boost::trim( _line );

            for (int _block = 0; _block < _n_blocks; _block++ ) {
                
                // first line gives the j index in the matrix
                //cout << _line << endl;
                
                boost::tokenizer<> tok( _line );
                std::transform( tok.begin(), tok.end(), std::back_inserter( _j_indeces ), 
                  &boost::lexical_cast<int,std::string> );
                //std::copy( _j_indeces.begin(), _j_indeces.end(), std::ostream_iterator<int>(std::cout,"\n") );
            
                // read the block of max _basis_size lines + the following header
                for (int i = 0; i <= _basis_size; i++ ) {
                    getline (_input_file, _line); 
                    //cout << _line << endl;
                    if ( std::string::npos == _line.find("D") ) break;
                    
                    // split the line on the i index and the rest
                    
                    vector<string> _row;
                    boost::trim( _line );
                    boost::algorithm::split( _row , _line, boost::is_any_of("\t "), boost::algorithm::token_compress_on); 
                   
                            
                    int _i_index = boost::lexical_cast<int>(  _row.front()  ); 
                    _row.erase( _row.begin() );
                    
                    //cout << _i_index << ":" << _line << endl ;
                    
                    std::vector<int>::iterator _j_iter = _j_indeces.begin();
                    
                    for (std::vector<string>::iterator iter = _row.begin()++; iter != _row.end(); iter++) {
                        string  _coefficient = *iter;
                       
                        boost::replace_first( _coefficient, "D", "e" );
                        //cout << boost::lexical_cast<double>( _coefficient ) << endl;
                        
                        int _j_index = *_j_iter;                                
                        _overlap( _i_index-1 , _j_index-1 ) = boost::lexical_cast<double>( _coefficient );
                        _j_iter++;
                        
                    }
 
                    
                }
                
                // clear the index for the next block
                _j_indeces.clear();        
            } // end of the blocks
            cout << "... finished reading the overlap matrix from " << filename << endl;
        } // end of the if "Overlap" found   
    } // end of the loop over the file
}


bool Orbitals::ParseGaussianLog( const char * filename ){
    
    string _line;
    unsigned _basis_size = 0;
    bool _read_overlap = false;
    
    cout << endl << "... parsing " << filename << endl;
    
    ifstream _input_file(filename);
    if (_input_file.fail()) {
        throw std::runtime_error("File filename is not found." );
        return 1;
    };
    
    while ( _input_file  ) {

        getline(_input_file, _line);
        //cout << _line << endl;
        
        // if a line has "electrons", must be the one we need
        std::string::size_type electrons_pos = _line.find("alpha electrons");
 
        if (electrons_pos != std::string::npos && !_has_electrons ) {          
            vector<string> results;
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t "),
                    boost::algorithm::token_compress_on); 
            
            _electrons = boost::lexical_cast<int>( results.front() );
            _has_electrons = true;
            cout << "... number of electrons: " << _electrons << endl;
 
        }
        
        std::string::size_type basis_pos = _line.find("basis functions");

        if (basis_pos != std::string::npos && !_has_basis_set_size ) {          
            vector<string> results;
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t "),
                    boost::algorithm::token_compress_on); 
            
            _basis_set_size = boost::lexical_cast<int>( results.front() );
            _has_basis_set_size = true;
            cout << "... basis set size: " << _basis_set_size << endl;
 
        }        
        
        // count occupied and unoccupied orbitals
        std::string::size_type eigenvalues_pos = _line.find("Alpha");

        while (eigenvalues_pos != std::string::npos && !_has_occupied_levels && !_has_unoccupied_levels) {    

            boost::trim( _line );            
            std::list<std::string> stringList;
            boost::iter_split(stringList, _line, boost::first_finder("--"));
            
            vector<string> energies; 
            boost::trim( stringList.back() );
            
            //cout << stringList.back() << endl;
                    
            boost::algorithm::split(energies, stringList.back(), boost::is_any_of("\t "),
            boost::algorithm::token_compress_on);            
            
            if  ( stringList.front().find("virt.") != std::string::npos ) {
                _unoccupied_levels += energies.size();
                //cout << energies.size() << endl;
                energies.clear();
            }
            
            if  ( stringList.front().find("occ.") != std::string::npos ) {
                _occupied_levels += energies.size();
                energies.clear();
            }
            
            getline(_input_file, _line);
            eigenvalues_pos = _line.find("Alpha");
            
            boost::iter_split(stringList, _line, boost::first_finder("--"));
            
            if ( eigenvalues_pos == std::string::npos ) {
                _has_occupied_levels = true;
                _has_unoccupied_levels = true;
                cout << "... ... occupied levels: " << _occupied_levels << endl;
                cout << "... ... unoccupied levels: " << _unoccupied_levels << endl;
            }
        }               
        // check if all information has been accumulated
        if ( _has_electrons && _has_basis_set_size && _has_occupied_levels ) break;
    }
    
    cout << "... done parsing " << filename << " file" << endl;
}

const int    &Orbitals::getBasisSetSize() const { 
    if ( _has_basis_set_size ) {
        return _basis_set_size; 
    } else {
        throw std::runtime_error(" Basis set size is unknown. Parse a log file first. " );
    }
}

const int     &Orbitals::getNumberOfLevels() const {
    if ( _has_occupied_levels && _has_unoccupied_levels ) {
        return  _occupied_levels + _unoccupied_levels; 
    } else {
        throw std::runtime_error(" Number of levels is unknown. Parse a log file first. " );
    }
}

const int     &Orbitals::getNumberOfElectrons() const {
    if ( _has_electrons ) {
        return  _electrons; 
    } else {
        throw std::runtime_error(" Number of electrons is unknown. Parse a log file first. " );
    }
}

/**
 * 
 * @param _energy_difference [ev] Two levels are degenerate if their energy is smaller than this value
 * @return A map with key as a level and a vector which is a list of close lying orbitals
 */
bool Orbitals::CheckDegeneracy( double _energy_difference ) {
    
    ub::vector<double>::iterator it1 = _mo_energies.begin();
    bool _degenerate = false;
    
    cout << endl <<"... checking level degeneracy " << endl;
    
    _level_degeneracy.clear();
            
    while ( it1 !=_mo_energies.end() ) {
        
        // in all containers counters start with 0; real life - with 1
        int _level1 = std::distance(_mo_energies.begin(), it1) + 1;
        
        // add the level itself - it is easier to loo over all levels later
        _level_degeneracy[_level1].push_back( _level1 );        
        
        ub::vector<double>::iterator it2 = it1;
        it2++;
        
        while (  it2 !=_mo_energies.end() ) {
            //cout << _level1 << ":" << *it1 << ":" << *it2 << endl;
            double energy1 = *it1;
            double energy2 = *it2;
            
            // in all containers counters start with 0; real life - with 1
            int _level2 = std::distance(_mo_energies.begin(), it2) + 1;
            
            if ( abs(energy1 - energy2)*_conv_Hrt_eV < _energy_difference ) {
                _level_degeneracy[_level1].push_back( _level2 );
                _level_degeneracy[_level2].push_back( _level1 );
                _degenerate = true;
            }
            it2++;
        }
        it1++;
    }

    if ( _degenerate ) {
        cout << "... ... some levels are degenerate" << endl; 
        for (std::map<int, std::vector<int> >::iterator it = _level_degeneracy.begin();
                it != _level_degeneracy.end();
                ++it) {
            // output only degenerate levels
            if ( (it->second).size() > 1 ) {
                std::cout << "....level  "<< it->first << " : ";
                std::vector<int> level_list;
                for (vector<int>::iterator lev = (it->second).begin(); lev != (it->second).end(); lev++)
                        cout << *lev << " ";
                cout << endl;
            }
        }
    } else {
        cout << "... ... no degeneracy found" << endl;  
    }
     
    cout << "... done checking level degeneracy" << endl;   
    
    _has_degeneracy = true;
    return _degenerate;
    
}    

std::vector<int>* Orbitals::getDegeneracy( int level, double _energy_difference ) {
    if ( !_has_degeneracy ) {
        
        CheckDegeneracy( _energy_difference );       
        /* 

        int _ld = _level_degeneracy.at(level).size();       
        if ( _ld > 1 ) {
                cout << "....level " << level << " degeneracy is: " <<  _ld << endl;
        } else {
                cout << "....level " << level << " is not degenerate" << endl;
        }
        */
        
    }        
    
    return &_level_degeneracy.at(level);
}

 bool Orbitals::Save( const char * filename ) {
        // create and open an archive for input
        //std::ifstream ifs( filename );
        //boost::archive::text_iarchive ia(ifs);
        //boost::serialization::guid_defined(*this);
        //ia << (*this);
 }

/*
 template<typename Archive> 
void Orbitals::serialize(Archive& ar, const unsigned version) {

    ar & _has_basis_set_size;
    ar & _has_occupied_levels;
    ar & _has_unoccupied_levels;
    ar & _has_electrons;
    ar & _has_degeneracy;
    
    //ar & _basis_set_size;
    //ar & _occupied_levels;
    //ar & _unoccupied_levels;
    //ar & _electrons;
       
    if ( _has_degeneracy ) {
        //ar & _level_degeneracy;
    }
    
    //std::vector<int>                    _active_levels;
    //ub::vector<double>                  _mo_energies;    
    //ub::matrix<double>                  _mo_coefficients;
    //ub::symmetric_matrix<double>        _overlap;
 
}*/
        
}}
