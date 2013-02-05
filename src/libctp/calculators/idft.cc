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
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
// +++++++++++++++++++++++++++++ //
// IDFT MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IDFT::Initialize(ctp::Topology *top, tools::Property* options ) {
    
    ParseOptionsXML( options );
    _orbitalsA.ReadOrbitalsGaussian( _orbitalsA_file.c_str() );
    _orbitalsAB.ReadOverlapGaussian( _logAB_file.c_str() );
    //SQRTOverlap();
    CalculateJ();
}

    
void IDFT::ParseOptionsXML( tools::Property *opt ) {
   
    // Orbitals are in fort.7 file; number of electrons in .log file
    
    // Molecule A
    string key = "options.idft.moleculeA";

    cout << key + ".orbitals" << endl;
    
    if ( opt->exists( key + ".orbitals" ) ) {
        _orbitalsA_file = opt->get( key + ".orbitals" ).as< string > ();
    }
    else {
        cout << key + ".orbitals" << endl;
        exit(0);
        throw std::runtime_error("Error in options: molecule A orbitals filename is missing.");
    }
    
    if ( opt->exists( key + ".log" ) ) {
        _logA_file = opt->get( key + ".log" ).as< string > ();
    }
    else {
        cout << key + ".orbitals" << endl;
        exit(0);
        throw std::runtime_error("Error in options: molecule A log filename is missing.");
    }
    
    // Molecule B
    key = "options.idft.moleculeB";

    if ( opt->exists( key+".orbitals" ) ) {
        _orbitalsB_file = opt->get( key + ".orbitals" ).as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: molecule B orbitals filename is missing.");
    }

    if ( opt->exists( key + ".log" ) ) {
        _logB_file = opt->get( key + ".log" ).as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: molecule B log filename is missing.");
    }
    
    // Dimer 
    key = "options.idft.moleculeAB";
    
    if ( opt->exists(key + ".orbitals") ) {
        _orbitalsAB_file = opt->get(key + ".orbitals").as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: dimer orbitals filename is missing.");
    }

    if ( opt->exists( key + ".log") ) {
        _logAB_file = opt->get(key + ".log").as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: dimer log filename is missing.");
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
class _inv_sqrt {
public:
  double operator()(double x) { return 1./x; }
};
*/

double inv_sqrt(double x) { return 1./sqrt(x); }

/*
 * Calculates S^{-1/2}
 */
void IDFT::SQRTOverlap() {
       
    double (*_inv_sqrt)(double);
    _inv_sqrt = &inv_sqrt;

    ub::vector<double>                  _eigenvalues;
    ub::matrix<double>                  _eigenvectors;
    ub::symmetric_matrix<double>        _overlap;

    _overlap = *_orbitalsAB.getOverlap();
    int _basis_size = _overlap.size1(); 

    cout << "Calculating SQRT of the " << _basis_size << "x" << _basis_size  << " overlap matrix" << endl;

    _eigenvalues.resize( _basis_size );
    _eigenvectors.resize( _basis_size, _basis_size ); 
    
    
//  test case  

/*
    int _basis_size = 3;
    _overlap.resize( _basis_size ); 
    _eigenvalues.resize( _basis_size );
    _eigenvectors.resize( _basis_size, _basis_size ); 
    
    //eigenvalues 3, 6, 9
    //eigenvectors (1,2,2), (-2,-1,2), (2,-2,1)
   
    _overlap(0,0) = 7;   
    _overlap(1,0) =-2;  _overlap(1,1) = 6;  
    _overlap(2,0) = 0;  _overlap(2,1) =-2; _overlap(2,2) = 5;

*/
    
    EigenvaluesSymmetric(_overlap, _eigenvalues, _eigenvectors);
    cout << "..eigenvalue problem solved " << endl;
    //cout << _eigenvalues << endl;
    //cout << _eigenvectors << endl;
     
    // compute inverse sqrt of all eigenvalues
    std::transform(_eigenvalues.begin(), _eigenvalues.end(), _eigenvalues.begin(),  _inv_sqrt );

    // form a diagonal matrix S^{-1/2}
    ub::diagonal_matrix<double> _diagS2( _eigenvalues.size(), _eigenvalues.data() ); 

    // multiply from the left on the U
    ub::matrix<double> _temp = ub::prod( _eigenvectors, _diagS2 );
    
    // multiply from the right on the transpose U
    ub::trans(_eigenvectors);
    ub::matrix<double> S2 = ub::prod( _temp, _eigenvectors);
    cout << "..projection matrix constructed  " << endl;
    
    
    
    /* for the test case above S2 has the following form 
    * [[0.3937418627,0.07087375404,0.0209304492],
    *  [0.07087375404,0.4501091889,0.0918042032],
    *  [0.0209304492,0.0918042032,0.4750808413]]
    */

    // cleanup
    _diagS2.clear();
    _temp.clear();

    //cout << S2 << endl;
    //cout << _overlap << endl;
    
    cout << "Done with the overlap matrix" << endl;
    
    
 }

/*
void IDFT::CleanUp() {

}
*/



void IDFT::CalculateJ() {

    // building the outer product of the dimer   
    ub::matrix<double> _monomersAB (4, 4);
    ub::matrix_range< ub::matrix<double> > mr (_monomersAB, ub::range (0, 2), ub::range (0, 2));
    for (unsigned i = 0; i < mr.size1 (); ++ i)
        for (unsigned j = 0; j < mr.size2 (); ++ j)
            mr (i, j) = 1;
    
    ub::matrix_range< ub::matrix<double> > mr1 (_monomersAB, ub::range (2, 4), ub::range (2, 4));
    for (unsigned i = 0; i < mr1.size1 (); ++ i)
        for (unsigned j = 0; j < mr1.size2 (); ++ j)
            mr1 (i, j) = 2;

    std::cout << _monomersAB << std::endl;
    
    ub::matrix<double> C(2, 2);
    C(0,0) = 3; C(0,1) = 3;
    C(1,0) = 3; C(1,1) = 3;
    project(_monomersAB, ub::range (2, 4), ub::range (2, 4)) = C;
 
    std::cout << _monomersAB << std::endl;
        
       // what I really need is this:  project(A, r1, r2) = C; 
    

}

}};