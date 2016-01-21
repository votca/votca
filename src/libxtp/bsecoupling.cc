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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/xtp/votca_xtp_config.h>

#include <votca/xtp/bsecoupling.h>
#include <votca/tools/linalg.h>

#include <boost/format.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include <boost/progress.hpp>

using boost::format;

namespace votca { namespace xtp {

namespace ub = boost::numeric::ublas;


/*
 * Calculates S^{-1/2}
 */
void BSECoupling::SQRTOverlap(ub::symmetric_matrix<double> &S, 
                          ub::matrix<double> &S2 ) {
       
    //double (*_inv_sqrt)(double);
    //_inv_sqrt = &inv_sqrt;

    ub::vector<double>                  _eigenvalues;
    ub::matrix<double>                  _eigenvectors;

    int _size = S.size1(); 

    _eigenvalues.resize( _size );
    _eigenvectors.resize( _size, _size ); 
    
    votca::tools::linalg_eigenvalues_symmetric(S, _eigenvalues, _eigenvectors);
    
    // compute inverse sqrt of all eigenvalues
    // std::transform(_eigenvalues.begin(), _eigenvalues.end(), _eigenvalues.begin(),  _inv_sqrt );

    // form a diagonal matrix S^{-1/2}
    ub::diagonal_matrix<double> _diagS2( _eigenvalues.size(), _eigenvalues.data() ); 

    // multiply from the left on the U
    ub::matrix<double> _temp = ub::prod( _eigenvectors, _diagS2 );
    
    // multiply from the right on the transpose U
    S2 = ub::prod( _temp, ub::trans( _eigenvectors ) );
    
 }

float BSECoupling::getSingletCouplingElement( int levelA, int levelB,  Orbitals* _orbitalsA,
    Orbitals* _orbitalsB, ub::matrix<float>* _JAB, double  _energy_difference ) {

    const float _conv_Ryd_eV = 27.21138386/2.0;
    int _levelsA = _orbitalsA->BSESingletEnergies().size();
    if ( _levelsA == 0  ){
      _levelsA = _orbitalsA->BSETripletEnergies().size();
    }
    return _JAB->at_element( levelA  , levelB + _levelsA ) * _conv_Ryd_eV;


    /*int _statesA = _orbitalsA->BSE;
    int _statesB = _orbitalsB->getNumberOfLevels();    
    
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
        return _JAB->at_element( levelA  , levelB + _levelsA ) * _conv_Ryd_eV;
    }*/
    // the  matrix should be symmetric, could also return this element
    // _JAB.at_element( _levelsA + levelB - 1  , levelA - 1 );
}

float BSECoupling::getTripletCouplingElement( int levelA, int levelB,  Orbitals* _orbitalsA,
    Orbitals* _orbitalsB, ub::matrix<float>* _JAB, double  _energy_difference ) {

    const float _conv_Ryd_eV = 27.21138386/2.0;
    int _levelsA = _orbitalsA->BSETripletEnergies().size();
    return _JAB->at_element( levelA  , levelB + _levelsA ) * _conv_Ryd_eV;


    /*int _statesA = _orbitalsA->BSE;
    int _statesB = _orbitalsB->getNumberOfLevels();    
    
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
        return _JAB->at_element( levelA  , levelB + _levelsA ) * _conv_Ryd_eV;
    }*/
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
 bool BSECoupling::CalculateCouplings_OLD(Orbitals* _orbitalsA, Orbitals* _orbitalsB, 
    Orbitals* _orbitalsAB, ub::matrix<float>* _JAB) {

    LOG(logDEBUG,*_pLog) << "Calculating exciton couplings" << flush;
    


        
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        LOG(logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
        return false;
    }

    // number of levels stored in monomers
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    boost::timer t; // start timing
    double _st = t.elapsed();
    
    
    
        
    // get exciton information of molecule A
    ub::matrix<float>& _bseA = _orbitalsA->BSESingletCoefficients();
    int _bseA_cmax = _orbitalsA->getBSEcmax();
    int _bseA_cmin = _orbitalsA->getBSEcmin();
    int _bseA_vmax = _orbitalsA->getBSEvmax();
    int _bseA_vmin = _orbitalsA->getBSEvmin();
    int _bseA_vtotal = _bseA_vmax - _bseA_vmin +1 ;
    int _bseA_ctotal = _bseA_cmax - _bseA_cmin +1 ;
    int _bseA_size   = _bseA_vtotal * _bseA_ctotal;
    int _bseA_exc    = _bseA.size2();
    
    LOG(logDEBUG,*_pLog)  << " molecule A has " << _bseA_exc << " excitons with dimension " << _bseA_size << flush;
    // now, two storage assignment matrices for two-particle functions
    ub::matrix<int> _combA;
    _combA.resize(_bseA_size,2);
    int _cnt = 0;
    for ( int _v = 0; _v < _bseA_vtotal; _v++){
        for ( int _c = 0; _c < _bseA_ctotal; _c++){
            _combA(_cnt,0) = _v;
            _combA(_cnt,1) = _bseA_vtotal + _c;
            _cnt++;
        }
    }
    
    
    // get exciton information of molecule B
    ub::matrix<float>& _bseB = _orbitalsB->BSESingletCoefficients();
    int _bseB_cmax = _orbitalsB->getBSEcmax();
    int _bseB_cmin = _orbitalsB->getBSEcmin();
    int _bseB_vmax = _orbitalsB->getBSEvmax();
    int _bseB_vmin = _orbitalsB->getBSEvmin();
    int _bseB_vtotal = _bseB_vmax - _bseB_vmin +1 ;
    int _bseB_ctotal = _bseB_cmax - _bseB_cmin +1 ;
    int _bseB_size   = _bseB_vtotal * _bseB_ctotal;
    int _bseB_exc    = _bseB.size2();
    
    LOG(logDEBUG,*_pLog)  << " molecule B has " << _bseB_exc << " excitons with dimension " << _bseB_size << flush;
    // now, two storage assignment matrices for two-particle functions
    ub::matrix<int> _combB;
    _combB.resize(_bseB_size,2);
    _cnt = 0;
    for ( int _v = 0; _v < _bseB_vtotal; _v++){
        for ( int _c = 0; _c < _bseB_ctotal; _c++){
            _combB(_cnt,0) = _bseA_vtotal + _bseA_ctotal + _v;
            _combB(_cnt,1) = _bseA_vtotal + _bseA_ctotal + _bseB_vtotal + _c;
            _cnt++;
        }
    }
    
    // get exciton information of pair AB
    ub::matrix<float>& _bseAB = _orbitalsAB->BSESingletCoefficients();
    int _bseAB_cmax = _orbitalsAB->getBSEcmax();
    int _bseAB_cmin = _orbitalsAB->getBSEcmin();
    int _bseAB_vmax = _orbitalsAB->getBSEvmax();
    int _bseAB_vmin = _orbitalsAB->getBSEvmin();
    int _bseAB_vtotal = _bseAB_vmax - _bseAB_vmin +1 ;
    int _bseAB_ctotal = _bseAB_cmax - _bseAB_cmin +1 ;
    int _bseAB_size   = _bseAB_vtotal * _bseAB_ctotal;
    int _bseAB_exc    = _bseAB.size2();
    
    LOG(logDEBUG,*_pLog)  << " dimer AB has " << _bseAB_exc << " excitons with dimension " << _bseAB_size << flush;
    // now, two storage assignment matrices for two-particle functions
    ub::matrix<int> _combAB;
    _combAB.resize(_bseAB_size,2);
    _cnt = 0;
    for ( int _v = 0; _v < _bseAB_vtotal; _v++){
        for ( int _c = 0; _c < _bseAB_ctotal; _c++){
            _combAB(_cnt,0) = _v;
            _combAB(_cnt,1) = _bseAB_vtotal + _c;
            _cnt++;
        }
    }
    
    
    LOG(logDEBUG,*_pLog) << "Levels stored: Basis A[" << _levelsA << ":" << _basisA << "]"
                                     << " B[" << _levelsB << ":" << _basisB << "]" << flush;
    
    // DFT levels can be reduced to those used in BSE
    _levelsA = _bseA_vtotal + _bseA_ctotal;
    _levelsB = _bseB_vtotal + _bseB_ctotal;
    LOG(logDEBUG,*_pLog) << " Levels used in BSE of molA: " << _bseA_vmin << " to " << _bseA_cmax << " total: " << _bseA_vtotal + _bseA_ctotal <<  flush;
    LOG(logDEBUG,*_pLog) << " Levels used in BSE of molB: " << _bseB_vmin << " to " << _bseB_cmax << " total: " << _bseB_vtotal + _bseB_ctotal <<  flush;
    
    
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
            << _psi_AxB.size2() << "]";    
    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    
     
    //ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = *_orbitalsA->getOrbitals();
    // ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = *_orbitalsB->getOrbitals(); 

    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = ub::project( *_orbitalsA->getOrbitals() , ub::range(_bseA_vmin, _bseA_cmax+1) , ub::range ( 0, _basisA ));
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = ub::project( *_orbitalsB->getOrbitals(), ub::range(_bseB_vmin, _bseB_cmax+1) , ub::range ( 0, _basisB )); 
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

    // row: monomer orbitals, col: dimer orbitals
    ub::matrix<double> _proj_check = ub::zero_matrix<double>(_levelsA + _levelsB, _psi_AxB_dimer_basis.size2() );
    
    for ( int i = 0; i < _levelsA + _levelsB; i++ ){
        // occupied dimer levels
        for ( int j = 0; j < _bseAB_vtotal ; j++ ){
            
            _proj_check(i,0) += _psi_AxB_dimer_basis(i,j)*_psi_AxB_dimer_basis(i,j);
            
        }
        
        // empty dimer levels
        for ( unsigned j = _bseAB_vtotal ; j < _psi_AxB_dimer_basis.size2() ; j++ ){
            
            _proj_check(i,1) += _psi_AxB_dimer_basis(i,j)*_psi_AxB_dimer_basis(i,j);
            
        }
 
        
        cout << " Monomer level projection (occ:virt:total)" << i << " : " << _proj_check(i,0) << " : " << _proj_check(i,1) << " : " << _proj_check(i,0) + _proj_check(i,1) << endl; 
        
        
    }
    
    
 
    
    
    
    
    
    
    
    
    
    // some more convenient storage
    
    ub::matrix<float> _kap;
    _kap.resize(_bseA_size,_bseAB_size);
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) );
            
        }
    }

    
    ub::matrix<float> _kbp;
    _kbp.resize(_bseB_size,_bseAB_size);
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) );
        }
    }

    
     LOG(logDEBUG,*_pLog)  << " projection of product functions(" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();    
    
     
     // now project the exciton functions (ugly!)
     // ub::matrix<float> _proj_excA;
     // _proj_excA.resize(_bseA_exc,_bseAB_exc);
     ub::matrix<float> _temp = ub::prod( _kap, _bseAB );
     ub::matrix<float> _proj_excA = ub::prod( ub::trans(_bseA), _temp);
     
     
     for ( unsigned i = 0 ; i < _proj_excA.size1() ; i++ ){
     
        float check = 0.0;
        
        for ( unsigned j = 0; j < _proj_excA.size2(); j++ ){
            
            check += _proj_excA(i,j) * _proj_excA(i,j) ;
        }
        
        cout << " Monomer exciton " << i << " norm " << check << endl;
        
     }
     
     
     
     _temp = ub::prod( _kbp, _bseAB );
     ub::matrix<float> _proj_excB = ub::prod( ub::trans(_bseB ), _temp);
     
     
     cout << "exciton projectors: " << _proj_excA.size1() << " : " << _proj_excA.size2();
     
     
     
     
     
     
     LOG(logDEBUG,*_pLog)  << " projection of exciton wave functions (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();  
     
     
     ub::matrix<float> _Hamiltonian_AB = ub::zero_matrix<float>(_bseAB_size,_bseAB_size );
     std::vector<float>& _bseAB_energies = _orbitalsAB->BSESingletEnergies();

     for ( int _i=0; _i<_bseAB_size; _i++){
         _Hamiltonian_AB(_i,_i) = _bseAB_energies[_i];
     }
     
 
     
     ub::matrix<float> _J_dimer = ub::zero_matrix<float>( _bseA_exc + _bseB_exc , _bseA_exc + _bseB_exc );
     ub::matrix<float> _S_dimer = ub::zero_matrix<float>( _bseA_exc + _bseB_exc , _bseA_exc + _bseB_exc );
     
     // setup J
     _temp = ub::prod( _Hamiltonian_AB , ub::trans(_proj_excA) );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excA, _temp ); // J_AA = proj_excA x H x trans(proj_excA)
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excB, _temp ); // J_BA = proj_excB x H x trans(proj_excA)
     
     _temp = ub::prod( _Hamiltonian_AB , ub::trans(_proj_excB) );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc )  ) = ub::prod( _proj_excA, _temp ); // J_AB = proj_excA x H x trans(proj_excB)
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc)  ) = ub::prod( _proj_excB, _temp ); // J_BB = proj_excB x H x trans(proj_excB)
     
     LOG(logDEBUG,*_pLog)  << " setup J (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();  

     
     // setup S
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excA, ub::trans(_proj_excA) ); // S_AA = proj_excA x trans(proj_excA)
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excB, ub::trans(_proj_excA)); // S_BA = proj_excB x trans(proj_excA)
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc )  ) = ub::prod( _proj_excA, ub::trans(_proj_excB)); // J_AB = proj_excA x H x trans(proj_excB)
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc)  ) = ub::prod( _proj_excB, ub::trans(_proj_excB) ); // J_BB = proj_excB x H x trans(proj_excB)
     
     LOG(logDEBUG,*_pLog)  << " setup S (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();  
     
     ub::vector<float> _S_eigenvalues; 
     linalg_eigenvalues( _S_eigenvalues, _S_dimer);

     ub::matrix<float> _diagS = ub::zero_matrix<float>( _bseA_exc + _bseB_exc  , _bseA_exc + _bseB_exc );
     for ( int _i =0; _i < _bseA_exc + _bseB_exc ; _i++){
         _diagS(_i,_i) = 1.0/sqrt(_S_eigenvalues[_i]);
     }
     _temp=ub::prod( _diagS, ub::trans(_S_dimer) ) ;
     ub::matrix<float> _transform = ub::prod( _S_dimer, _temp );
     LOG(logDEBUG,*_pLog)  << " transformation matrix (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();  
     // final coupling elements
     _temp=ub::prod(_J_dimer, _transform);
     ub::matrix<float> _J_eff = ub::prod( _transform, _temp);
     LOG(logDEBUG,*_pLog)  << " effective couplings (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();
     
     
     LOG(logDEBUG,*_pLog)  << " singlet coupling: " << _J_eff(0,_bseA_exc+1)*13.6058 << " and " <<  _J_eff(_bseA_exc+1, 0) * 13.6058 << endl;
     LOG(logDEBUG,*_pLog)  << " singlet coupling: " << _J_eff(0,_bseA_exc)*13.6058 << " and " <<  _J_eff(_bseA_exc, 0) * 13.6058 << endl; 
     LOG(logDEBUG,*_pLog)  << " singlet coupling: " << _J_eff(1,_bseA_exc+1)*13.6058 << " and " <<  _J_eff(_bseA_exc+1, 1) * 13.6058 << endl; 


  
    LOG(logDEBUG,*_pLog)  << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed(); 
    
    LOG(logDEBUG,*_pLog) << "Done with exciton couplings" << flush;
    return true;   

} 

/**
 * \brief evaluates electronic couplings  
 * 
 * This is a different version with not diagonalized BSE Hamiltonian
 * 
 * @param _orbitalsA molecular orbitals of molecule A
 * @param _orbitalsB molecular orbitals of molecule B
 * @param _orbitalsAB molecular orbitals of the dimer AB
 * @param _JAB matrix with electronic couplings
 * @return false if failed
 */
bool BSECoupling::CalculateCouplings(Orbitals* _orbitalsA, Orbitals* _orbitalsB, 
    Orbitals* _orbitalsAB, ub::matrix<float>* _JAB_singlet, ub::matrix<float>* _JAB_triplet, string _type) {
          
    LOG(logDEBUG,*_pLog) << "  Calculating exciton couplings" << flush;
        

        
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        LOG(logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
        return false;
    }

    // number of levels stored in monomers
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    boost::timer t; // start timing
    //double _st = t.elapsed();
        
    // get exciton information of molecule A
    int _bseA_cmax        = _orbitalsA->getBSEcmax();
    int _bseA_cmin        = _orbitalsA->getBSEcmin();
    int _bseA_vmax        = _orbitalsA->getBSEvmax();
    int _bseA_vmin        = _orbitalsA->getBSEvmin();
    int _bseA_vtotal      = _bseA_vmax - _bseA_vmin +1 ;
    int _bseA_ctotal      = _bseA_cmax - _bseA_cmin +1 ;
    int _bseA_size        = _bseA_vtotal * _bseA_ctotal;
    int _bseA_singlet_exc = _orbitalsA->BSESingletCoefficients().size2();
    int _bseA_triplet_exc = _orbitalsA->BSETripletCoefficients().size2();

    LOG(logDEBUG,*_pLog)  << "   molecule A has " << _bseA_singlet_exc << " singlet excitons with dimension " << _bseA_size << flush;
    LOG(logDEBUG,*_pLog)  << "   molecule A has " << _bseA_triplet_exc << " triplet excitons with dimension " << _bseA_size << flush;
    
    // now, two storage assignment matrices for two-particle functions
    ub::matrix<int> _combA;
    _combA.resize(_bseA_size,2);
    int _cnt = 0;
    for ( int _v = 0; _v < _bseA_vtotal; _v++){
        for ( int _c = 0; _c < _bseA_ctotal; _c++){
            _combA(_cnt,0) = _v;
            _combA(_cnt,1) = _bseA_vtotal + _c;
            _cnt++;
        }
    }
    
    
    // get exciton information of molecule B
    int _bseB_cmax        = _orbitalsB->getBSEcmax();
    int _bseB_cmin        = _orbitalsB->getBSEcmin();
    int _bseB_vmax        = _orbitalsB->getBSEvmax();
    int _bseB_vmin        = _orbitalsB->getBSEvmin();
    int _bseB_vtotal      = _bseB_vmax - _bseB_vmin +1 ;
    int _bseB_ctotal      = _bseB_cmax - _bseB_cmin +1 ;
    int _bseB_size        = _bseB_vtotal * _bseB_ctotal;
    int _bseB_singlet_exc = _orbitalsB->BSESingletCoefficients().size2();
    int _bseB_triplet_exc = _orbitalsB->BSETripletCoefficients().size2();
    
    LOG(logDEBUG,*_pLog)  << "   molecule B has " << _bseB_singlet_exc << " singlet excitons with dimension " << _bseB_size << flush;
    LOG(logDEBUG,*_pLog)  << "   molecule B has " << _bseB_triplet_exc << " triplet excitons with dimension " << _bseB_size << flush;
    
    // now, two storage assignment matrices for two-particle functions
    ub::matrix<int> _combB;
    _combB.resize(_bseB_size,2);
    _cnt = 0;
    for ( int _v = 0; _v < _bseB_vtotal; _v++){
        for ( int _c = 0; _c < _bseB_ctotal; _c++){
            _combB(_cnt,0) = _bseA_vtotal + _bseA_ctotal + _v;
            _combB(_cnt,1) = _bseA_vtotal + _bseA_ctotal + _bseB_vtotal + _c;
            _cnt++;
        }
    }
    
    // get exciton information of pair AB
    int _bseAB_cmax = _orbitalsAB->getBSEcmax();
    int _bseAB_cmin = _orbitalsAB->getBSEcmin();
    int _bseAB_vmax = _orbitalsAB->getBSEvmax();
    int _bseAB_vmin = _orbitalsAB->getBSEvmin();
    int _bseAB_vtotal = _bseAB_vmax - _bseAB_vmin +1 ;
    int _bseAB_ctotal = _bseAB_cmax - _bseAB_cmin +1 ;
    int _bseAB_size   = _bseAB_vtotal * _bseAB_ctotal;
    // check if electron-hole interaction matrices are stored
    if ( ! _orbitalsAB->hasEHinteraction() ){
        LOG(logERROR,*_pLog) << "BSE EH int not stored in dimer " << flush;
        return false;
    }
    ub::matrix<float>&    _eh_d = _orbitalsAB->eh_d(); 
    ub::matrix<float>&    _eh_x = _orbitalsAB->eh_x(); 
    
    LOG(logDEBUG,*_pLog)  << "   dimer AB has BSE EH interaction (direct)   with dimension " << _eh_d.size1() << " x " <<  _eh_d.size2() << flush;
    LOG(logDEBUG,*_pLog)  << "   dimer AB has BSE EH interaction (exchange) with dimension " << _eh_x.size1() << " x " <<  _eh_x.size2() << flush;
    // now, two storage assignment matrices for two-particle functions
    ub::matrix<int> _combAB;
    _combAB.resize(_bseAB_size,2);
    _cnt = 0;
    for ( int _v = 0; _v < _bseAB_vtotal; _v++){
        for ( int _c = 0; _c < _bseAB_ctotal; _c++){
            //_combAB(_cnt,0) = _v;
            //_combAB(_cnt,1) = _bseAB_vtotal + _c;
            
            _combAB(_cnt,0) = _bseAB_vmin + _v;
            _combAB(_cnt,1) = _bseAB_vmin + _bseAB_vtotal + _c;
            
            _cnt++;
        }
    }
    
    // DFT levels can be reduced to those used in BSE
    _levelsA = _bseA_vtotal + _bseA_ctotal;
    _levelsB = _bseB_vtotal + _bseB_ctotal;
    LOG(logDEBUG,*_pLog) << "   levels used in BSE of molA: " << _bseA_vmin << " to " << _bseA_cmax << " total: " << _bseA_vtotal + _bseA_ctotal <<  flush;
    LOG(logDEBUG,*_pLog) << "   levels used in BSE of molB: " << _bseB_vmin << " to " << _bseB_cmax << " total: " << _bseB_vtotal + _bseB_ctotal <<  flush;
    
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        LOG(logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
        return false;
    } 
    
    //       | Orbitals_A          0 |      | Overlap_A |     
    //       | 0          Orbitals_B |  X   | Overlap_B |  X  Transpose( Orbitals_AB )
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
    ub::matrix<double> _psi_AxB ( _levelsA + _levelsB, _basisA + _basisB  );
    

    // constructing merged orbitals
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = ub::project( *_orbitalsA->getOrbitals() , ub::range(_bseA_vmin, _bseA_cmax+1) , ub::range ( 0, _basisA ));
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = ub::project( *_orbitalsB->getOrbitals(), ub::range(_bseB_vmin, _bseB_cmax+1) , ub::range ( 0, _basisB )); 
    
    // psi_AxB * S_AB * psi_AB
    LOG(logDEBUG,*_pLog) << "   projecting monomer onto dimer orbitals" << flush; 
    ub::matrix<double> _orbitalsAB_Transposed = ub::trans( *_orbitalsAB->getOrbitals() );  
    if ( (*_orbitalsAB->getOverlap()).size1() == 0 ) {
            LOG(logERROR,*_pLog) << "Overlap matrix is not stored"; 
            return false;
    }
   
    ub::matrix<double> _psi_AB = ub::prod( *_orbitalsAB->getOverlap(), _orbitalsAB_Transposed );  
    ub::matrix<double> _psi_AxB_dimer_basis = ub::prod( _psi_AxB, _psi_AB );  
    _psi_AB.clear();

    // _psi_AxB_dimer_basis = T in notes, dimension ( LA + LB, LD)
    
    
    //cout << "Size of _psi_AxB_dimer_basis " << _psi_AxB_dimer_basis.size1() << " : " <<  _psi_AxB_dimer_basis.size2() << flush; 
    
      
    // some more convenient storage
    ub::matrix<float> _kap;
    _kap.resize(_bseA_size,_bseAB_size);
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) );
            
        }
    }

    

    
    
    ub::matrix<float> _kbp;
    _kbp.resize(_bseB_size,_bseAB_size);
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) );
        }
    }
    
    
  

    LOG(logDEBUG,*_pLog)  << "   construct projection of product functions " << flush; 

    LOG(logDEBUG,*_pLog)  << "   Evaluating " << _type << flush; 
    
        
    //     cout << "Size of _kap " << _kap.size1() << " : " <<  _kap.size2() << "\n" << flush; 
    //     cout << "Size of _kbp " << _kbp.size1() << " : " <<  _kbp.size2() << "\n" << flush; 
    
    // now the different spin types
    if ( _type == "singlets" || _type == "all"){
      // get singlet BSE Hamiltonian from _orbitalsAB
      ub::matrix<float> _Hamiltonian_AB = _eh_d + 2.0 * _eh_x;
        ub::matrix<float>& _bseA = _orbitalsA->BSESingletCoefficients();
        ub::matrix<float>& _bseB = _orbitalsB->BSESingletCoefficients();
        
      //  cout << "Size of _Hamiltonian_AB " << _Hamiltonian_AB.size1() << " : " <<  _Hamiltonian_AB.size2() << flush;
// cout << "Size of _bseA " << _bseA.size1() << " : " <<  _bseA.size2() << flush; 
   //      cout << "Size of _bseB " << _bseB.size1() << " : " <<  _bseB.size2() << flush; 

     
        bool _singlets = ProjectExcitons( _kap, _kbp, _bseA, _bseB, _Hamiltonian_AB, *_JAB_singlet);
        if ( _singlets ) {
            LOG(logDEBUG,*_pLog)  << "   calculated singlet couplings " << flush;
        }
        /*
        // diagonalize the effective Hamiltonian?
        ub::vector<float> _coupled_energies;
        ub::matrix<float> _coupled_coefficients;
	std::vector< std::vector<double> >& _free_dipolesA = _orbitalsA->TransitionDipoles();
	std::vector< std::vector<double> >& _free_dipolesB = _orbitalsB->TransitionDipoles();
        linalg_eigenvalues(*_JAB_singlet, _coupled_energies, _coupled_coefficients,_JAB_singlet->size1() );
        LOG(logDEBUG,*_pLog)  << "   calculated EVs of coupling matrix " << flush;
	cout << "\n" << endl;
        for ( int i =0 ; i< _JAB_singlet->size1(); i++){
            
            
	    std::vector<double> tdipole(3,0.0);
            for ( int j = 0; j < 50; j++){
	      tdipole[0] += _coupled_coefficients(j,i)*_free_dipolesA[j][0] + _coupled_coefficients(j+50,i)*_free_dipolesB[j][0];
	      tdipole[1] += _coupled_coefficients(j,i)*_free_dipolesA[j][1] + _coupled_coefficients(j+50,i)*_free_dipolesB[j][1];
	      tdipole[2] += _coupled_coefficients(j,i)*_free_dipolesA[j][2] + _coupled_coefficients(j+50,i)*_free_dipolesB[j][2];
	    }
	    double tdipole_strength = tdipole[0]*tdipole[0] + tdipole[1]*tdipole[1] + tdipole[2]*tdipole[2];
            double oscillator_strength = tdipole_strength * _coupled_energies(i) /3.0;

	    LOG(logINFO, *_pLog) << (format("  S = %1$4d Omega = %2$+1.4f eV  lamdba = %3$+3.2f nm ") % (i + 1) % (13.6058 * _coupled_energies(i)) % (1240.0/(13.6058 * _coupled_energies(i))) ).str() << flush;
	    LOG(logINFO, *_pLog) << (format("           TrDipole length gauge   dx = %1$+1.4f dy = %2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f") % (tdipole[0]) % (tdipole[1]) % (tdipole[2]) % (tdipole_strength) % (oscillator_strength)).str() << flush;
	    for (int _i_bse = 0; _i_bse < _JAB_singlet->size1(); _i_bse++) {
	      // if contribution is larger than 0.2, print
	      double _weight = pow(_coupled_coefficients(_i_bse, i), 2);
	      if (_weight > 0.2) {
                if ( _i_bse < 50) {
		LOG(logINFO, *_pLog) << (format("           EXCITON A %1$-3d : %2$3.1f%%") % _i_bse % (100.0 * _weight)).str() << flush;
		} else
		  {
		    LOG(logINFO, *_pLog) << (format("           EXCITON B %1$-3d : %2$3.1f%%") % (_i_bse-50) % (100.0 * _weight)).str() << flush;
		  }
	      }
	    }
	    LOG(logINFO, *_pLog) << (format("   ")).str() << flush;




            cout << " E" << i << " : " << _coupled_energies(i)*13.605 << " TD " << tdipole[0] << " " << tdipole[1] << " " << tdipole[2] << endl;          
	}
        //LOG(logDEBUG,*_pLog)  << " singlet coupling: " << _JAB_singlet->at_element(0,_bseA_singlet_exc)*13.6058 << " and " <<  _JAB_singlet->at_element(_bseA_singlet_exc, 0) * 13.6058 << endl; 
	*/
    }
    
                
        
    if ( _type == "triplets" || _type == "all"){
        LOG(logDEBUG,*_pLog)  << "   Evaluating " << _type << flush; 
        // get triplet BSE Hamiltonian from _orbitalsAB
        ub::matrix<float> _Hamiltonian_AB = _eh_d;
        
     /*   cout << endl;
        for ( int i = 0; i<_eh_d.size1(); i++ ){
            for ( int j = 0; j<_eh_d.size2(); j++ ){
                
                
                cout << "eh D [" << i << " : " << j << "]: " << _eh_d(i,j) << endl; 
                
            }
            
        }*/
        
        
        
        ub::matrix<float>& _bseA = _orbitalsA->BSETripletCoefficients();
        ub::matrix<float>& _bseB = _orbitalsB->BSETripletCoefficients();
        
      /*   cout << endl;
        for ( int i = 0; i<_bseA.size1(); i++ ){
            for ( int j = 0; j<_bseA.size2(); j++ ){
                
                
                cout << "BSE A and B [" << i << " : " << j << "]: " << _bseA(i,j) << " and " << _bseB(i,j)  << endl; 
                
            }
        
        
        }*/
        
        
        bool _triplets = ProjectExcitons( _kap, _kbp, _bseA, _bseB, _Hamiltonian_AB, *_JAB_triplet);
        if ( _triplets ) {
            LOG(logDEBUG,*_pLog)  << "   calculated triplet couplings " << flush;
        }
        LOG(logDEBUG,*_pLog)  << " triplet coupling: " << _JAB_triplet->at_element(0,_bseA_triplet_exc)*13.6058 << " and " <<  _JAB_triplet->at_element(_bseA_triplet_exc, 0) * 13.6058 << endl; 
    }
    
    
    LOG(logDEBUG,*_pLog) << "  Done with exciton couplings" << flush;
    return true;   
};


bool BSECoupling::ProjectExcitons(ub::matrix<float>& _kap, ub::matrix<float>& _kbp, ub::matrix<float>& _bseA, ub::matrix<float>& _bseB, ub::matrix<float>& _H, ub::matrix<float>& _J){
    
    
    cout << " Dimensions of _bseA " << _bseA.size1() << " : " << _bseA.size2() << endl;
    
     // get projection of monomer excitons on dimer product functions
     ub::matrix<float> _proj_excA = ub::prod( ub::trans( _bseA ), _kap);
     ub::matrix<float> _proj_excB = ub::prod( ub::trans( _kbp ), _bseB);

      
     int _bseA_exc = _bseA.size2();
     int _bseB_exc = _bseB.size2();
     ub::matrix<float> _J_dimer = ub::zero_matrix<float>( _bseA_exc + _bseB_exc , _bseA_exc + _bseB_exc );
     ub::matrix<float> _S_dimer = ub::zero_matrix<float>( _bseA_exc + _bseB_exc , _bseA_exc + _bseB_exc );
     
     // setup J
     ub::matrix<float> _temp = ub::prod( _H , ub::trans(_proj_excA) );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excA, _temp ); // J_AA = proj_excA x H x trans(proj_excA)
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( ub::trans(_proj_excB), _temp ); // J_BA = proj_excB x H x trans(proj_excA)
     
     _temp = ub::prod( _H , _proj_excB );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc )  ) = ub::prod( _proj_excA, _temp ); // J_AB = proj_excA x H x trans(proj_excB)
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc)  ) = ub::prod( ub::trans(_proj_excB), _temp ); // J_BB = proj_excB x H x trans(proj_excB)


     cout << "J_dimer " << _J_dimer(0,50)*13.605 << " and " << _J_dimer(50,0)*13.605 << endl;
     
     
     // setup S
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excA, ub::trans(_proj_excA) ); // S_AA = proj_excA x trans(proj_excA)
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( 0, _bseA_exc )  ) = ub::trans( ub::prod( _proj_excA, _proj_excB)); // S_BA = proj_excB x trans(proj_excA)
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc )  ) = ub::prod( _proj_excA, _proj_excB); // S_AB = proj_excA x H x trans(proj_excB)
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc)  ) = ub::prod( ub::trans(_proj_excB),_proj_excB ); // S_BB = proj_excB x H x trans(proj_excB)

     ub::vector<float> _S_eigenvalues; 

//    ub::matrix<double> _S_dimer_double = _S_dimer;
  //  ub::vector<double> _S_eigenvalues_double; 
     
     
     linalg_eigenvalues( _S_eigenvalues, _S_dimer);
    // linalg_eigenvalues( _S_eigenvalues_double, _S_dimer_double);

     //cout << endl;
     //for ( int i = 0 ; i < _S_eigenvalues.size(); i++){
     //    cout << " S eigenvalues " << _S_eigenvalues[i] << " and " <<  _S_eigenvalues_double[i] << endl; 
     // }
     
     

     

     if ( _S_eigenvalues[0] < 0.0 ) {
         
         cerr << " \n Negative eigenvalues in projected overlap. Projection quality not sufficient. Try increasing dimer basis. " << endl;
         return false;
         
     }

     
     
     ub::matrix<float> _diagS = ub::zero_matrix<float>( _bseA_exc + _bseB_exc  , _bseA_exc + _bseB_exc );
     for ( int _i =0; _i < _bseA_exc + _bseB_exc ; _i++){

         _diagS(_i,_i) = 1.0/sqrt(_S_eigenvalues[_i]);
     }

 
    
     
     
     
     
     //ub::matrix<float> _transform = ub::prod( _S_dimer, ub::prod( _diagS, ub::trans(_S_dimer) )  );

     ub::matrix<float> _transtemp = ub::prod( _diagS, ub::trans(_S_dimer));
     ub::matrix<float> _transform = ub::prod( _S_dimer,_transtemp );

     ub::matrix<float> _J_temp = ub::prod(_J_dimer, _transform);
     
             // final coupling elements
     // _J = ub::prod( _transform, ub::prod(_J_dimer, _transform));
    _J = ub::prod( _transform, _J_temp);
    
         //cout << endl;
     for ( unsigned i = 0 ; i < _J.size1(); i++){
         for ( unsigned j = 0 ; j < _J.size1(); j++){
           cout << " J [" << i << " : " << j << "] " << _J(i,j)*13.605 << " and " <<  _J_dimer(i,j)*13.605 << endl; 
         }
         }  
    
    
    
    
     return true;
    
    
    
    
}
    
}}
