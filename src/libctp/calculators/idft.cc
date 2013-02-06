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
    _orbitalsB.ReadOrbitalsGaussian( _orbitalsB_file.c_str() );
    _orbitalsAB.ReadOrbitalsGaussian( _orbitalsAB_file.c_str() );
    
    _orbitalsAB.ReadOverlapGaussian( _logAB_file.c_str() );
    
    //SQRTOverlap();
    
    _orbitalsA.ParseGaussianLog(_logA_file.c_str());
    _orbitalsB.ParseGaussianLog(_logB_file.c_str());
    _orbitalsAB.ParseGaussianLog(_logAB_file.c_str());    
    
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

    cout << "S2: " << S2 << endl;
    cout << "Overlap: " << _overlap << endl;
    
    cout << "Done with the overlap matrix" << endl;
    
    
 }

/*
void IDFT::CleanUp() {

}
*/



void IDFT::CalculateJ() {
 
    double conv_Hrt_eV=27.21138386;
            
    /* test case
    ub::matrix<double> _monomersAB (4, 5);
    ub::zero_matrix<double> _AB (4, 5);

    _monomersAB = _AB;
    
    std::cout << _monomersAB << std::endl;
    
    ub::matrix<double> C(2, 2);
    C(0,0) = 3; C(0,1) = 3;
    C(1,0) = 3; C(1,1) = 3;
    
    ub::matrix<double> B(2, 2);
    B(0,0) = 5; B(0,1) = 5;
    B(1,0) = 5; B(1,1) = 5;
    
    ub::project(_monomersAB, ub::range (2, 4), ub::range (3, 5)) = C;
    ub::project(_monomersAB, ub::range (0, 2), ub::range (0, 2)) = B;

    std::cout << _monomersAB << std::endl;
    */
    
    
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA.getBasisSetSize();
    int _basisB = _orbitalsB.getBasisSetSize();
    
    cout << "basis [A:B] " << _basisA << ":" << _basisB << endl;
    
    int _levelsA = _orbitalsA.getNumberOfLevels();
    int _levelsB = _orbitalsB.getNumberOfLevels();
    
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
        
    //cout << zeroB << endl;
    //cout << zeroA << endl;
    
    ub::matrix<double> _psi_AxB ( _levelsA + _levelsB, _basisA + _basisB  );
    
    // AxB = | A 0 |  //
    //       | 0 B |  //      
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;
    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = *_orbitalsA.getOrbitals();
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = *_orbitalsB.getOrbitals();    
    //cout << "_psi_AxB: " << _psi_AxB << endl;
    
    // Fock matrix of a dimer   
    ub::diagonal_matrix<double> _fock_AB( _orbitalsAB.getNumberOfLevels(), (*_orbitalsAB.getEnergies()).data() ); 
    //cout << "_fock_AxB: " << _psi_AxB << endl;
    
    //cout << "Overlap: "  << *_orbitalsAB.getOverlap();
    //cout << "OrbitalsAB: "  << *_orbitalsAB.getOrbitals();
    
    // psi_AxB * SAB * psi_AB
    ub::matrix<double> _psi_AB = ub::prod( *_orbitalsAB.getOverlap(), ub::trans( *_orbitalsAB.getOrbitals() ) );   
         
    ub::matrix<double> _psi_AxB_dimer_basis = ub::prod( _psi_AxB, _psi_AB );

    //cout << (*_orbitalsAB.getOverlap()).at_element(0,0) << endl;
    //cout << (*_orbitalsAB.getOrbitals()).at_element(0,0) << endl;
    //cout << _psi_AB.at_element(0,0) << endl;
    
    /*
    for (int i = 0; i < _psi_AxB_dimer_basis.size1(); i++ ) {
        for (int j = 0; j < _psi_AxB_dimer_basis.size2(); j++ ) {
            cout << i << " " << j << " " << _psi_AxB_dimer_basis.at_element(i, j) << endl;
            
        }
    }
    exit(0);
     */
    
    //cout << "_psi_AxB_dimer_basis: " << _psi_AxB_dimer_basis << endl;
    _psi_AB.clear();
    
    // J = psi_AxB_dimer_basis * FAB * psi_AxB_dimer_basis^T
    ub::matrix<double> _temp = ub::prod( _fock_AB, ub::trans( _psi_AxB_dimer_basis ) ) ;
    //cout << "_temp: " << _temp << endl;
    ub::matrix<double> JAB = ub::prod( _psi_AxB_dimer_basis, _temp);
    _temp.clear();
    _fock_AB.clear();
    
    // S = psi_AxB_dimer_basis * psi_AxB_dimer_basis^T
    ub::matrix<double> _S = ub::prod( _psi_AxB_dimer_basis, ub::trans( _psi_AxB_dimer_basis ));
    cout << "SAxB: " << _S << endl;
    
    
    //cout << "JAB: " << JAB << endl;
    int HOMO_A = _orbitalsA.getNumberOfElectrons() - 1 ;
    int HOMO_B = _orbitalsB.getNumberOfElectrons() - 1 ;
    
    cout << JAB.at_element( HOMO_A , HOMO_B + _levelsA ) * conv_Hrt_eV << endl; 
    cout << JAB.at_element(_levelsA + HOMO_B, HOMO_A ) * conv_Hrt_eV << endl;
    
}

}};