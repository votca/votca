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
};   

void Orbitals::Initialize( tools::Property *options )
{
    //this->ParseOrbitalsXML( options );
}

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
        cerr << "File " << filename << " with molecular orbitals is not found " << endl;
        return 1;
    };

    // number of coefficients per line is  in the first line of the file (5D15.8)
    getline(_input_file, _line);
    std::vector<string> strs;
    boost::algorithm::split(strs, _line, boost::is_any_of("(D)"));
    //cout << strs.at(1) << endl;
    int nrecords_in_line = boost::lexical_cast<int>(strs.at(1));
    string format = strs.at(2);

    cout << "Orbital file: "
            << nrecords_in_line << " records per line, in D"
            << format << " format." << endl;

    while (_input_file) {

        getline(_input_file, _line);
        // if a line has an equality sign, must be energy
        std::string::size_type energy_pos = _line.find("=");

        if (energy_pos != std::string::npos) {

            vector<string> results;
            boost::algorithm::split(results, _line, boost::is_any_of("\t ="),
                    boost::algorithm::token_compress_on);
            //cout << results[1] << ":" << results[2] << ":" << results[3] << ":" << results[4] << ":" << results[5] << endl;

            _level = boost::lexical_cast<int>(results[1]);
            boost::replace_first(results[5], "D", "e");
            _energies[ _level ] = boost::lexical_cast<double>(results[5]);            
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
   
   return 0;
}




}}
