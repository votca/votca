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


#include "idft.h"
#include "votca/ctp/qmcalculator.h"
#include <votca/ctp/eigenvalues.h>
#include <boost/numeric/ublas/io.hpp>

namespace votca { namespace ctp {
    
// +++++++++++++++++++++++++++++ //
// IDFT MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IDFT::Initialize(ctp::Topology *top, tools::Property* options ) {
    
    ParseOptionsXML( options );
    //_orbitalsA.ReadOrbitalsGaussian( _orbitalsA_file.c_str() );
    //_orbitalsAB.ReadOverlapGaussian( _overlapAB_file.c_str() );
    SQRTOverlap();
}

    
void IDFT::ParseOptionsXML( tools::Property *opt ) {

    string key = "options.idft";
    
    if ( opt->exists(key+".orbitals_A") ) {
        _orbitalsA_file = opt->get(key + ".orbitals_A").as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: molecule A orbitals filename is missing.");
    }
    
    if ( opt->exists(key+".orbitals_B") ) {
        _orbitalsB_file = opt->get(key + ".orbitals_B").as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: molecule B orbitals filename is missing.");
    }
   
    if ( opt->exists(key+".orbitals_AB") ) {
        _overlapAB_file = opt->get(key + ".overlap_AB").as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: dimer orbitals filename is missing.");
    }
   
    /* --- ORBITALS.XML Structure ---
     * <options>
     *   <idft>
     *     <orbitals_A>fort.7</orbitals_A>
     *     <orbitals_B>fort.7</orbitals_B>
     *     <orbitals_AB>fort.7</orbitals_AB>
     *     <overlap_AB>dimer.log</overlap_AB>
     *   </idft>
     * </options>
     */

}

/*
 * Calculates S^{-1/2}
 */
void IDFT::SQRTOverlap() {
    
    boost::numeric::ublas::vector<double>       _eigenvalues;
    boost::numeric::ublas::matrix<double>       _eigenvectors;
    boost::numeric::ublas::matrix<double>       _overlap;
 
    //_overlap = _orbitalsAB._overlap;
    
    int _basis_size = _orbitalsAB.getBasisSetSize(); 
    
    _overlap.resize( _basis_size, _basis_size ); 
    _eigenvalues.resize( _basis_size );
    _eigenvectors.resize( _basis_size, _basis_size ); 
    
/* test case 
    eigenvalues 3, 6, 9
    eigenvectors (1,2,2), (-2,-1,2), (2,-2,1)

    _overlap(0,0) = 7;  _overlap(0,1) =-2; _overlap(0,2) = 0; 
    _overlap(1,0) =-2;  _overlap(1,1) = 6; _overlap(1,2) =-2; 
    _overlap(2,0) = 0;  _overlap(2,1) =-2; _overlap(2,2) = 5;
    
*/
    
    EigenvaluesSymmetric(_overlap, _eigenvalues, _eigenvectors);
    
    cout << _eigenvalues << endl;
    cout << _eigenvectors << endl;
}

/*
void IDFT::CleanUp() {

}
*/


/*
void IDFT::CalculateJ(QMPair *pair) {

}
*/
    
    
}};