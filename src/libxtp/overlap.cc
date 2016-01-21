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

#define NDEBUG
#define OVERLAP_DEBUG

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/xtp/votca_xtp_config.h>

#include <votca/xtp/overlap.h>
#include <votca/tools/linalg.h>

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include <boost/progress.hpp>

namespace votca { namespace xtp {

namespace ub = boost::numeric::ublas;

double inv_sqrt(double x) { return 1./sqrt(x); }

/*
 * Calculates S^{-1/2}
 */
void Overlap::SQRTOverlap(ub::symmetric_matrix<double> &S, 
                          ub::matrix<double> &S2 ) {
       
    double (*_inv_sqrt)(double);
    _inv_sqrt = &inv_sqrt;

    ub::vector<double>                  _eigenvalues;
    ub::matrix<double>                  _eigenvectors;

    int _size = S.size1(); 

    _eigenvalues.resize( _size );
    _eigenvectors.resize( _size, _size ); 
    
    votca::tools::linalg_eigenvalues_symmetric(S, _eigenvalues, _eigenvectors);
    
    // compute inverse sqrt of all eigenvalues
    std::transform(_eigenvalues.begin(), _eigenvalues.end(), _eigenvalues.begin(),  _inv_sqrt );

    // form a diagonal matrix S^{-1/2}
    ub::diagonal_matrix<double> _diagS2( _eigenvalues.size(), _eigenvalues.data() ); 

    // multiply from the left on the U
    ub::matrix<double> _temp = ub::prod( _eigenvectors, _diagS2 );
    
    // multiply from the right on the transpose U
    S2 = ub::prod( _temp, ub::trans( _eigenvectors ) );
    
 }

double Overlap::getCouplingElement( int levelA, int levelB,  Orbitals* _orbitalsA,
    Orbitals* _orbitalsB, ub::matrix<double>* _JAB, double  _energy_difference ) {

    const double _conv_Hrt_eV = 27.21138386;   

            int _levelsA = _orbitalsA->getNumberOfLevels();
    //int _levelsB = _orbitalsB->getNumberOfLevels();    
    
    if ( _energy_difference != 0 ) {
        std::vector<int> list_levelsA = *_orbitalsA->getDegeneracy( levelA, _energy_difference );
        std::vector<int> list_levelsB = *_orbitalsA->getDegeneracy( levelB, _energy_difference );
        
        double _JAB_sq = 0; double _JAB_one_level;
        
        for (std::vector<int>::iterator iA = list_levelsA.begin()++; iA != list_levelsA.end(); iA++) {
                for (std::vector<int>::iterator iB = list_levelsB.begin()++; iB != list_levelsB.end(); iB++) { 
                    //cout << *iA << ':' << *iB << endl;
                    _JAB_one_level = _JAB->at_element( *iA - 1  , *iB -1 + _levelsA );
                    _JAB_sq +=  _JAB_one_level*_JAB_one_level ;
                }
        }
        
        return sqrt(_JAB_sq / ( list_levelsA.size() * list_levelsB.size() ) ) * _conv_Hrt_eV ;
        
    } else {
        return _JAB->at_element( levelA - 1  , levelB -1 + _levelsA ) * _conv_Hrt_eV;
    }
    // the  matrix should be symmetric, could also return this element
    // _JAB.at_element( _levelsA + levelB - 1  , levelA - 1 );
}

/**
 * \brief evaluates electronic couplings  
 * 
 * This is a slow version with a rather large block matrix AxB
 * 
 * @param _orbitalsA molecular orbitals of molecule A
 * @param _orbitalsB molecular orbitals of molecule B
 * @param _orbitalsAB molecular orbitals of the dimer AB
 * @param _JAB matrix with electronic couplings
 * @return false if failed
 */
bool Overlap::CalculateIntegrals(Orbitals* _orbitalsA, Orbitals* _orbitalsB, 
    Orbitals* _orbitalsAB, ub::matrix<double>* _JAB) {

    LOG(logDEBUG,*_pLog) << "Calculating electronic couplings" << flush;
        
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        LOG(logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
        return false;
    }
        
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    boost::timer t; // start timing
    double _st = t.elapsed();
    
    LOG(logDEBUG,*_pLog) << "Levels:Basis A[" << _levelsA << ":" << _basisA << "]"
                                     << " B[" << _levelsB << ":" << _basisB << "]" << flush;
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        LOG(logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
        return false;
    } 
    
    //       | Orbitals_A          0 |      | Overlap_A |     
    //       | 0          Orbitals_B |  X   | Overlap_B |  X  Transpose( Orbitals_AB )
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
    ub::matrix<double> _psi_AxB ( _levelsA + _levelsB, _basisA + _basisB  );
    

     LOG(logDEBUG,*_pLog) << "Constructing direct product AxB [" 
            << _psi_AxB.size1() << "x" 
            << _psi_AxB.size2() << "]" ;    
    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = *_orbitalsA->getOrbitals();
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = *_orbitalsB->getOrbitals(); 

    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();

    
    // psi_AxB * S_AB * psi_AB
    LOG(logDEBUG,*_pLog) << "Projecting dimer onto monomer orbitals"; 
    ub::matrix<double> _orbitalsAB_Transposed = ub::trans( *_orbitalsAB->getOrbitals() );  
    if ( (*_orbitalsAB->getOverlap()).size1() == 0 ) {
            LOG(logERROR,*_pLog) << "Overlap matrix is not stored"; 
            return false;
    }
    #ifdef OVERLAP_DEBUG 
        cout << "\n\t\tprod1 [" 
             << (*_orbitalsAB->getOverlap()).size1() << "x" << (*_orbitalsAB->getOverlap()).size2() << "] ["
             << _orbitalsAB_Transposed.size1()  << "x" << _orbitalsAB_Transposed.size2() << "] ";  
    #endif    
    ub::matrix<double> _psi_AB = ub::prod( *_orbitalsAB->getOverlap(), _orbitalsAB_Transposed );  
    #ifdef OVERLAP_DEBUG 
        cout << "\t\tprod2 [" 
             << _psi_AxB.size1() << "x" << _psi_AxB.size2() << "] ["
             << _psi_AB.size1()  << "x" << _psi_AB.size2() << "] ";  
    #endif     
    ub::matrix<double> _psi_AxB_dimer_basis = ub::prod( _psi_AxB, _psi_AB );  
    _psi_AB.clear();
    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();    

     
    // J = psi_AxB_dimer_basis * FAB * psi_AxB_dimer_basis^T
    LOG(logDEBUG,*_pLog) << "Projecting the Fock matrix onto the dimer basis";   
    ub::diagonal_matrix<double> _fock_AB( _orbitalsAB->getNumberOfLevels(), (*_orbitalsAB->getEnergies()).data() ); 
    #ifdef OVERLAP_DEBUG 
        cout << "\n\t\tprod3 [" 
             << _fock_AB.size1() << "x" << _fock_AB.size2() << "] T["
             << _psi_AxB_dimer_basis.size1()  << "x" << _psi_AxB_dimer_basis.size2() << "] ";  
    #endif
    ub::matrix<double> _temp = ub::prod( _fock_AB, ub::trans( _psi_AxB_dimer_basis ) ) ; 
    #ifdef OVERLAP_DEBUG 
        cout << "\t\tprod4 [" 
             << _psi_AxB_dimer_basis.size1() << "x" << _psi_AxB_dimer_basis.size2() << "] ["
             << _temp.size1()  << "x" << _temp.size2() << "] ";  
    #endif   
    ub::matrix<double> JAB_dimer = ub::prod( _psi_AxB_dimer_basis, _temp);  
    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();    

    // S = psi_AxB_dimer_basis * psi_AxB_dimer_basis^T
    LOG(logDEBUG,*_pLog) << "Constructing Overlap matrix";    
    #ifdef OVERLAP_DEBUG 
        cout << "\n\t\tprod5 [" 
             << _psi_AxB_dimer_basis.size1() << "x" << _psi_AxB_dimer_basis.size2() << "] T["
             << _psi_AxB_dimer_basis.size1()  << "x" << _psi_AxB_dimer_basis.size2() << "] ";  
    #endif
    ub::symmetric_matrix<double> _S_AxB = ub::prod( _psi_AxB_dimer_basis, ub::trans( _psi_AxB_dimer_basis ));  
    ub::matrix<double> _S_AxB_2(_S_AxB.size1(), _S_AxB.size1() );
    ub::trans( _S_AxB );
    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();    

    // Square root of the overlap matrix
    LOG(logDEBUG,*_pLog) << "Calculating square root of the overlap matrix";    
    SQRTOverlap( _S_AxB , _S_AxB_2 );
    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();    

     
    LOG(logDEBUG,*_pLog) << "Calculating the effective overlap JAB [" 
              << JAB_dimer.size1() << "x" 
              << JAB_dimer.size2() << "]";  
       
    ub::matrix<double> JAB_temp( _levelsA + _levelsB, _levelsA + _levelsB ); 
    #ifdef OVERLAP_DEBUG 
        cout << "\n\t\tprod6 [" 
             << JAB_dimer.size1() << "x" << JAB_dimer.size2() << "] T["
             << _S_AxB_2.size1()  << "x" << _S_AxB_2.size2() << "] ";  
    #endif
    ub::noalias(JAB_temp) = ub::prod( JAB_dimer, _S_AxB_2 );  
    #ifdef OVERLAP_DEBUG
        cout << "\t\tprod7 [" 
             << _S_AxB_2.size1() << "x" << _S_AxB_2.size2() << "] T["
             << JAB_temp.size1()  << "x" << JAB_temp.size2() << "] ";  
    #endif
    (*_JAB) = ub::prod( _S_AxB_2, JAB_temp );    
    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed(); 
    
    LOG(logDEBUG,*_pLog) << "Done with electronic couplings" << flush;
    return true;   

}

/**
 * \brief evaluates electronic couplings  
 * 
 * This is a different version with block matrices [slower]
 * 
 * @param _orbitalsA molecular orbitals of molecule A
 * @param _orbitalsB molecular orbitals of molecule B
 * @param _orbitalsAB molecular orbitals of the dimer AB
 * @param _JAB matrix with electronic couplings
 * @return false if failed
 */
bool Overlap::CalculateIntegralsOptimized(Orbitals* _orbitalsA, Orbitals* _orbitalsB, 
    Orbitals* _orbitalsAB, ub::matrix<double>* _JAB) {
          
    LOG(logDEBUG,*_pLog) << "Calculating electronic couplings" << flush;
        
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        LOG(logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
        return false;
    }
        
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    boost::timer t; // start timing
    double _st = t.elapsed();
    
    LOG(logDEBUG,*_pLog) << "Levels:Basis A[" << _levelsA << ":" << _basisA << "]"
                                     << " B[" << _levelsB << ":" << _basisB << "]" << flush;
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        LOG(logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
        return false;
    } 
     
    // these flags should be set before any ublas header is called 
    // #define NDEBUG
    // otherwise the code is very inefficient
    
    //       | Orbitals_A          0 |      | Overlap_A |     
    //       | 0          Orbitals_B |  X   | Overlap_B |  X  Transpose( Orbitals_AB )
    //constructing a slice of the Overlap matrix
    ub::matrix_range< ub::symmetric_matrix<double> > Overlap_A = ub::project( *_orbitalsAB->getOverlap(), ub::range ( 0, _basisA), ub::range (0, _basisA +_basisB) );
    ub::matrix_range< ub::symmetric_matrix<double> > Overlap_B = ub::project( *_orbitalsAB->getOverlap(), ub::range ( _basisA, _basisA +_basisB ), ub::range (0, _basisA +_basisB) );
    
    LOG(logDEBUG,*_pLog) << "Projecting the monomer onto dimer orbitals [" << _levelsA + _levelsB << "x" << _basisA + _basisB << "]";   
    ub::matrix<double> _psi_AB ( _levelsA + _levelsB, _basisA + _basisB  );

    ub::matrix_range< ub::matrix<double> > _psi_AB_A = ub::project( _psi_AB, ub::range (0, _levelsA ), ub::range ( 0, _basisA +_basisB ) ) ;
    ub::noalias(_psi_AB_A) = ub::prod(*_orbitalsA->getOrbitals(), Overlap_A);

    ub::matrix_range< ub::matrix<double> > _psi_AB_B = ub::project( _psi_AB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA +_basisB ) ) ;
    ub::noalias(_psi_AB_B) = ub::prod(*_orbitalsB->getOrbitals(), Overlap_B );
    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();
 
    ub::matrix<double> _psi_AxB_dimer_basis (_levelsA + _levelsB, _basisA + _basisB );
    ub::matrix<double> OrbAB_Transp = ub::trans( *_orbitalsAB->getOrbitals() );
    LOG(logDEBUG,*_pLog)  << "Transposing OrbitalsAB (" << t.elapsed() - _st << "s)" << "\x1b[0;39m" << flush; _st = t.elapsed();
    ub::noalias(_psi_AxB_dimer_basis) = ub::prod( _psi_AB,  OrbAB_Transp );
  
    LOG(logDEBUG,*_pLog)  << "Multiplying PsiAB x OrbitalsAB (" << t.elapsed() - _st << "s)" << "\x1b[0;39m" << flush; _st = t.elapsed();
    
    _psi_AB.resize(0,0,false); OrbAB_Transp.resize(0,0,false);
    
    //   _psi_AxB_dimer_basis * F  * _psi_AxB_dimer_basis^T
    LOG(logDEBUG,*_pLog) << "Projecting the Fock matrix onto the dimer basis";    
    ub::zero_matrix<double> _zero ( _levelsA + _levelsB, _levelsA + _levelsB );
    ub::matrix<double> JAB_dimer( _zero ) ;
    ub::vector<double> energies = (*_orbitalsAB->getEnergies());

    for ( int i1 = 0; i1 < _levelsA + _levelsB ; i1++ ) {
    for ( int i2 = i1; i2 < _levelsA + _levelsB; i2++ ) {
        for ( int k = 0; k < _basisA + _basisB; k++  ) {
                JAB_dimer(i1,i2) += _psi_AxB_dimer_basis.at_element(i1, k) * _psi_AxB_dimer_basis.at_element(i2, k) * energies(k);
        }
        JAB_dimer(i2,i1) = JAB_dimer(i1,i2);
    }}   
    energies.clear();
    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s)" << "\x1b[0;39m" << flush; _st = t.elapsed();
    
    // S = psi_AxB_dimer_basis * psi_AxB_dimer_basis^T
    ub::symmetric_matrix<double> _S_AxB = ub::prod( _psi_AxB_dimer_basis, ub::trans( _psi_AxB_dimer_basis ));
    _psi_AxB_dimer_basis.resize(0,0,false);

    ub::matrix<double> _S_AxB_2(_S_AxB.size1(), _S_AxB.size1() );
    
     ub::trans( _S_AxB );
     LOG(logDEBUG,*_pLog) << "Calculating square root of the overlap matrix [" 
             << _S_AxB.size1() << "x" 
             << _S_AxB.size2() << "]";    
     SQRTOverlap( _S_AxB , _S_AxB_2 );        
     _S_AxB.resize(0,0,false); 
     LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s)" << "\x1b[0;39m" << flush; _st = t.elapsed();
    
    
     LOG(logDEBUG,*_pLog) << "Calculating the effective overlap JAB [" 
             << JAB_dimer.size1() << "x" 
             << JAB_dimer.size2() << "]";    
    
    ub::matrix<double> JAB_temp( _levelsA + _levelsB, _levelsA + _levelsB );
    
    ub::noalias(JAB_temp) = ub::prod( JAB_dimer, _S_AxB_2 );
    (*_JAB) = ub::prod( _S_AxB_2, JAB_temp );
    
    // cleanup
    JAB_dimer.resize(0,0,false); JAB_temp.resize(0,0,false); _S_AxB_2.resize(0,0,false);
    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s)" << "\x1b[0;39m" << flush; _st = t.elapsed();
    
    LOG(logDEBUG,*_pLog) << "Done with electronic couplings" << flush;
    return true;   
};

    
}}
