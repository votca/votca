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


namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
// +++++++++++++++++++++++++++++ //
// IDFT MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IDFT::Initialize(ctp::Topology *top, tools::Property* options ) {
    
    _has_integrals = false;
    _has_degeneracy = false;
    
    ParseOptionsXML( options );
    
    //_orbitalsA.ReadOrbitalsGaussian( _orbitalsA_file.c_str() );
    //_orbitalsB.ReadOrbitalsGaussian( _orbitalsB_file.c_str() );
    //_orbitalsAB.ReadOrbitalsGaussian( _orbitalsAB_file.c_str() );
    
    //_orbitalsA.ParseGaussianLog(_logA_file.c_str());
    //_orbitalsB.ParseGaussianLog(_logB_file.c_str());
    //_orbitalsAB.ParseGaussianLog(_logAB_file.c_str());    

    //_orbitalsAB.ReadOverlapGaussian( _logAB_file.c_str() );
    
    int HOMO_A = _orbitalsA.getNumberOfElectrons() ;
    int HOMO_B = _orbitalsB.getNumberOfElectrons() ;
    
    int LUMO_A = _orbitalsA.getNumberOfElectrons() + 1;
    int LUMO_B = _orbitalsB.getNumberOfElectrons() + 1;
    
    cout << getCouplingElement( HOMO_A , HOMO_B ) << endl; 
    cout << getCouplingElement( LUMO_A , LUMO_B ) << endl; 

}

    
void IDFT::ParseOptionsXML( tools::Property *opt ) {
   
    // Orbitals are in fort.7 file; number of electrons in .log file
    
    string key = "options.idft.";   
    if ( opt->exists( key + ".degeneracy" ) ) {
        _has_degeneracy = true;
        _energy_difference = opt->get( key + ".degeneracy" ).as< double > ();
    }
    else {
        cout << "....not treating degenerate orbitals" << endl ;
    }    
    
    // Molecule A
    key = "options.idft.moleculeA.orbitals";

    //cout << key + ".orbitals" << endl;
    
    if ( opt->exists( key + ".file" ) ) {
        _orbitalsA_file = opt->get( key + ".file" ).as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: molecule A orbitals filename is missing.");
    }
    
    if ( opt->exists( key + ".occupied" ) ) {
        _max_occupied_levels = opt->get( key + ".occupied" ).as< int > ();
    }
    else {
        throw std::runtime_error("Maximum of occupied levels for molecule A is not provided");
    }
    
    key = "options.idft.moleculeA";

    if ( opt->exists( key + ".log" ) ) {
        _logA_file = opt->get( key + ".log" ).as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: molecule A log filename is missing.");
    }
    
    // Molecule B
    key = "options.idft.moleculeB.orbitals";

    if ( opt->exists( key+".file" ) ) {
        _orbitalsB_file = opt->get( key + ".file" ).as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: molecule B orbitals filename is missing.");
    }


    if ( opt->exists( key + ".occupied" ) ) {
        _max_occupied_levels = opt->get( key + ".occupied" ).as< int > ();
    }
    else {
        throw std::runtime_error("Maximum of occupied levels for molecule B is not provided");
    }

    key = "options.idft.moleculeB.";
    
    if ( opt->exists( key + ".log" ) ) {
        _logB_file = opt->get( key + ".log" ).as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: molecule B log filename is missing.");
    }
    
    // Dimer 
    key = "options.idft.moleculeAB.orbitals";
    
    if ( opt->exists(key + ".file") ) {
        _orbitalsAB_file = opt->get(key + ".file").as< string > ();
    }
    else {
        throw std::runtime_error("Error in options: dimer orbitals filename is missing.");
    }
    
    key = "options.idft.moleculeAB";
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
void IDFT::SQRTOverlap(ub::symmetric_matrix<double> &S, ub::matrix<double> &S2 ) {
       
    double (*_inv_sqrt)(double);
    _inv_sqrt = &inv_sqrt;

    ub::vector<double>                  _eigenvalues;
    ub::matrix<double>                  _eigenvectors;

    int _size = S.size1(); 

    cout << "....calculating SQRT of the " << _size << "x" << _size  << " overlap matrix" << endl;

    _eigenvalues.resize( _size );
    _eigenvectors.resize( _size, _size ); 
    
    
//  test case  

/*  test of the eigendecomposition code
    int _basis_size = 3;
    ub::symmetric_matrix<double> _overlap(_basis_size);
    _overlap.resize( _basis_size ); 
    _eigenvalues.resize( _basis_size );
    _eigenvectors.resize( _basis_size, _basis_size ); 
    
    //eigenvalues 3, 6, 9
    //eigenvectors (1,2,2), (-2,-1,2), (2,-2,1)
   
    _overlap(0,0) = 7;   
    _overlap(1,0) =-2;  _overlap(1,1) = 6;  
    _overlap(2,0) = 0;  _overlap(2,1) =-2; _overlap(2,2) = 5;

    EigenvaluesSymmetric(_overlap, _eigenvalues, _eigenvectors);
    cout << "....eigenvalue problem solved " << endl;
    
    cout << "eigenvalues" << _eigenvalues << endl;
    cout << "eigenvectors" << _eigenvectors << endl;
    
    ub::diagonal_matrix<double> _diag( _eigenvalues.size(), _eigenvalues.data() );
    ub::matrix<double> _left = ub::prod( _eigenvectors, _diag );
    cout <<  ub::prod( _left, ub::trans( _eigenvectors ) );
    
    exit(0);
*/    
    /* for the test case above S2 has the following form 
    * [[0.3937418627,0.07087375404,0.0209304492],
    *  [0.07087375404,0.4501091889,0.0918042032],
    *  [0.0209304492,0.0918042032,0.4750808413]]
    */
    
    
    EigenvaluesSymmetric(S, _eigenvalues, _eigenvectors);
    cout << "....eigenvalue problem solved " << endl;
    
    //cout << "eigenvalues" << _eigenvalues << endl;
    //cout << _eigenvectors << endl;     
    
    // compute inverse sqrt of all eigenvalues
    std::transform(_eigenvalues.begin(), _eigenvalues.end(), _eigenvalues.begin(),  _inv_sqrt );

    // form a diagonal matrix S^{-1/2}
    ub::diagonal_matrix<double> _diagS2( _eigenvalues.size(), _eigenvalues.data() ); 

    // multiply from the left on the U
    ub::matrix<double> _temp = ub::prod( _eigenvectors, _diagS2 );
    
    // multiply from the right on the transpose U
    S2 = ub::prod( _temp, ub::trans( _eigenvectors ) );
    cout << "....projection matrix constructed  " << endl;
       


    // cleanup
    _diagS2.clear();
    _temp.clear();

    //cout << "S2: " << S2 << endl;
    //cout << "Overlap: " << _overlap << endl;
    
    cout << "....done with the sqrt of a matrix" << endl;
    
    
 }

void IDFT::CalculateIntegrals() {
            
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
    
    cout << endl << "..calculating electronic couplings " << endl;
    
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA.getBasisSetSize();
    int _basisB = _orbitalsB.getBasisSetSize();
    
    
    cout << "....basis [molA:molB] " << _basisA << ":" << _basisB << endl;
    
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

    // psi_AxB * S_AB * psi_AB
    ub::matrix<double> _psi_AB = ub::prod( *_orbitalsAB.getOverlap(), ub::trans( *_orbitalsAB.getOrbitals() ) );          
    ub::matrix<double> _psi_AxB_dimer_basis = ub::prod( _psi_AxB, _psi_AB );
    _psi_AB.clear();
   
    /*
    for (int i = 0; i < _psi_AxB_dimer_basis.size1(); i++ ) {
        for (int j = 0; j < _psi_AxB_dimer_basis.size2(); j++ ) {
            cout << i << " " << j << " " << _psi_AxB_dimer_basis.at_element(i, j) << endl;
            
        }
    }
    exit(0);
     */
    
    // J = psi_AxB_dimer_basis * FAB * psi_AxB_dimer_basis^T
    ub::matrix<double> _temp = ub::prod( _fock_AB, ub::trans( _psi_AxB_dimer_basis ) ) ;
    ub::matrix<double> JAB_dimer = ub::prod( _psi_AxB_dimer_basis, _temp);
    _temp.clear(); _fock_AB.clear();
    
    // S = psi_AxB_dimer_basis * psi_AxB_dimer_basis^T
    ub::symmetric_matrix<double> _S_AxB = ub::prod( _psi_AxB_dimer_basis, ub::trans( _psi_AxB_dimer_basis ));
    //cout << "SAxB: " << _S_AxB << endl;

    /* test of an assignment 
    ub::matrix<double> C(2,2);
    C(0,0) = 1; C(0,1) = 2;
    C(1,0) = 2; C(1,1) = 3;
    
    ub::symmetric_matrix<double> B = C;
    cout << C << endl; 
    cout << B << endl; 
    */
       
    ub::matrix<double> _S_AxB_2(_S_AxB.size1(), _S_AxB.size1() );
    
    /* test of the SQRT routine
    ub::symmetric_matrix<double> _test(3,3);
    ub::matrix<double> _test2(3,3); 
    _test(0,0) = 7;   
    _test(1,0) =-2;  _test(1,1) = 6;  
    _test(2,0) = 0;  _test(2,1) =-2; _test(2,2) = 5; 
    SQRTOverlap(_test, _test2 );
    cout << _test2;
    exit(0);
    */
    ub::trans( _S_AxB );
    SQRTOverlap( _S_AxB , _S_AxB_2 );        
    _S_AxB.clear(); 
     
    cout << "....calculating the effective overlap"<< endl;
    ub::matrix<double> JAB_temp = prod( JAB_dimer, _S_AxB_2 );
    _JAB = ub::prod( _S_AxB_2, JAB_temp );
    
    // cleanup
    JAB_dimer.clear(); JAB_temp.clear(); _S_AxB_2.clear();
    
    //cout << JAB << endl;
    
    //cout << _S_AxB << endl;
    _has_integrals = true;
    cout << "..done calculating electronic couplings"<< endl;
       
    //cout << JAB_dimer.at_element( HOMO_A , HOMO_B + _levelsA ) * conv_Hrt_eV << endl; 
    //cout << JAB_dimer.at_element(_levelsA + HOMO_B, HOMO_A ) * conv_Hrt_eV << endl;

}

double IDFT::getCouplingElement( int levelA, int levelB ) {
    
    int _levelsA = _orbitalsA.getNumberOfLevels();
    int _levelsB = _orbitalsB.getNumberOfLevels();
    
    if ( !_has_integrals ) {
         CalculateIntegrals();       
    }
    
    if ( _has_degeneracy ) {
        std::vector<int> list_levelsA = *_orbitalsA.getDegeneracy( levelA, _energy_difference );
        std::vector<int> list_levelsB = *_orbitalsA.getDegeneracy( levelB, _energy_difference );
        
        double _JAB_sq = 0; double _JAB_one_level;
        
        for (std::vector<int>::iterator iA = list_levelsA.begin()++; iA != list_levelsA.end(); iA++) {
                for (std::vector<int>::iterator iB = list_levelsB.begin()++; iB != list_levelsB.end(); iB++) { 
                    //cout << *iA << ':' << *iB << endl;
                    _JAB_one_level = _JAB.at_element( *iA - 1  , *iB -1 + _levelsA );
                    _JAB_sq +=  _JAB_one_level*_JAB_one_level ;
                }
        }
        
        return sqrt(_JAB_sq / ( list_levelsA.size() * list_levelsB.size() ) ) * _conv_Hrt_eV ;
        
    } else {
        return _JAB.at_element( levelA - 1  , levelB -1 + _levelsA ) * _conv_Hrt_eV;
    }
    // the  matrix should be symmetric, could also return this element
    // _JAB.at_element( _levelsA + levelB - 1  , levelA - 1 );
}
    
/*
void IDFT::CleanUp() {

}
*/

}};