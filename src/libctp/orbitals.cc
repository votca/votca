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
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>

using namespace std;

namespace votca { namespace ctp {

Orbitals::Orbitals() { 
     
};   
    
Orbitals::~Orbitals() { 
    _mo_energies.clear();
    _mo_coefficients.clear();
    _overlap.clear();
};   

 /*
 * Reads in the MO coefficients from a GAUSSIAN fcheck file
 */
bool Orbitals::ReadOrbitalsGaussian( const char * filename )
{
    map <int, vector<double> > _coefficients;
    map <int, double> _energies;
    
    string _line;
    unsigned _levels = 0;
    unsigned _level;
    unsigned _basis_size = 0;

    ifstream _input_file(filename);
    if (_input_file.fail()) {
        cerr << endl << "File " << filename << " with molecular orbitals is not found " << endl;
        return 1;
    };

    // number of coefficients per line is  in the first line of the file (5D15.8)
    getline(_input_file, _line);
    std::vector<string> strs;
    boost::algorithm::split(strs, _line, boost::is_any_of("(D)"));
    //cout << strs.at(1) << endl;
    int nrecords_in_line = boost::lexical_cast<int>(strs.at(1));
    string format = strs.at(2);

    cout << endl << "Orbital file " << filename << " has " 
            << nrecords_in_line << " records per line, in D"
            << format << " format." << endl;

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
    cout << "Read in " << _levels << " levels. ";

    std::map< int, vector<double> >::iterator iter = _coefficients.begin();
    _basis_size = iter->second.size();

    for (iter = _coefficients.begin()++; iter != _coefficients.end(); iter++) {
        if (iter->second.size() != _basis_size) {
            cerr << "Error reading " << filename << ". Basis set size change from level to level.";

        }
    }
    cout << "Basis set size: " << _basis_size << "." << endl;

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
   
   cout << "Done reading orbital files from " << filename << endl;
   return 0;
}


 /*
 * Reads in the Orbital Overlap matrix from a GAUSSIAN log file
 */
bool Orbitals::ReadOverlapGaussian( const char * filename )
{
    string _line;
    unsigned _basis_size = 0;
    bool _read_overlap = false;
    
    cout << "Reading the overlap matrix from " << filename << endl;
    
    ifstream _input_file(filename);
    if (_input_file.fail()) {
        cerr << endl << "File " << filename << " with overlap matrix is not found " << endl;
        return 1;
    };
    
    while (_input_file && !_read_overlap ) {

        getline(_input_file, _line);
        // if a line has an equality sign, must be energy
        std::string::size_type energy_pos = _line.find("Overlap");
        std::string::size_type nbasis_pos = _line.find("NBasis");
 
        if (nbasis_pos != std::string::npos && _basis_size == 0 ) {
          
            vector<string> results;
            boost::trim( _line );
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                    boost::algorithm::token_compress_on); 
            //cout << results[1] << ":" << results[2] << ":" << results[3] << ":" << results[4] << endl;
            
            _basis_size = boost::lexical_cast<int>(results[1]);
            cout << "Number of basis functions: " << _basis_size << endl;
            // preparing the matrix container
            _overlap.resize( _basis_size );  
 
        }

        if (energy_pos != std::string::npos ) {
            
            _read_overlap = true;
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
        } // end of the if "Overlap" found       
    } // end of the loop over the file
    //cout << _overlap << endl;
   
}

}}
