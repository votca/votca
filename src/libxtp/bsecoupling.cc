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
#include <votca/tools/constants.h>
#include <boost/format.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include <boost/progress.hpp>



namespace votca { namespace xtp {

namespace ub = boost::numeric::ublas;
using boost::format;

void BSECoupling::Initialize(Property* options){
    
    std::string key = "options." + Identify(); 
    _doSinglets=false;
    _doTriplets=false;
    
    
    _spintype   = options->get(key + ".type").as<string> ();
        if(_spintype=="all"){
            _doSinglets=true;
            _doTriplets=true;
        }
        else if(_spintype=="triplets"){
            _doTriplets=true;
        }
        else if(_spintype=="singlets"){
            _doSinglets=true;
        }
        else{
            throw std::runtime_error((boost::format("Choice % for type not known.") % _spintype).str());
        }
    _degeneracy = options->get(key + ".degeneracy").as<double> ();
    
        _levA  = options->get(key + ".moleculeA.states").as<int> ();
        _levB  = options->get(key + ".moleculeB.states").as<int> ();
        _occA  = options->get(key + ".moleculeA.occLevels").as<int> ();
        _occB  = options->get(key + ".moleculeB.occLevels").as<int> ();
        _unoccA  = options->get(key + ".moleculeA.unoccLevels").as<int> ();
        _unoccB  = options->get(key + ".moleculeB.unoccLevels").as<int> ();
        
}

void BSECoupling::addoutput(Property *_type_summary,Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB){
    if (_doSinglets){
        Property *_singlet_summary = &_type_summary->add("singlets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB <_levB ; ++stateB ) {
               float JAB = getSingletCouplingElement( stateA , stateB );
               float energyAD = getSingletDimerEnergy( stateA  );
               float energyBD = getSingletDimerEnergy( stateB  );
               Property *_coupling_summary = &_singlet_summary->add("coupling", boost::lexical_cast<string>(JAB)); 
               float energyA = _orbitalsA->BSESingletEnergies()[stateA]*conv::ryd2ev_f;
               float energyB = _orbitalsB->BSESingletEnergies()[stateB]*conv::ryd2ev_f;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", energyA);
               _coupling_summary->setAttribute("energyB", energyB);
               _coupling_summary->setAttribute("energyA_dimer", energyAD);
               _coupling_summary->setAttribute("energyB_dimer", energyBD);
           } 
        }
    }
    if ( _doTriplets){
        Property *_triplet_summary = &_type_summary->add("triplets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB < _levA ; ++stateB ) {
               float JAB = getTripletCouplingElement( stateA , stateB );
               float energyAD = getTripletDimerEnergy( stateA  );
               float energyBD = getTripletDimerEnergy( stateB  );
               Property *_coupling_summary = &_triplet_summary->add("coupling", boost::lexical_cast<string>(JAB)); 
               float energyA = _orbitalsA->BSETripletEnergies()[stateA]*conv::ryd2ev_f;
               float energyB = _orbitalsB->BSETripletEnergies()[stateB]*conv::ryd2ev_f;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", energyA);
               _coupling_summary->setAttribute("energyB", energyB);
               _coupling_summary->setAttribute("energyA_dimer", energyAD);
               _coupling_summary->setAttribute("energyB_dimer", energyBD);
           } 
        }
    }       
}



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

float BSECoupling::getSingletCouplingElement( int levelA, int levelB) {

    

    return JAB_singlet( levelA  , levelB +  _levA ) * votca::tools::conv::ryd2ev_f;
}
float BSECoupling::getSingletDimerEnergy( int level) {

    

    return JAB_singlet( level  , level ) * votca::tools::conv::ryd2ev_f;
}
float BSECoupling::getTripletDimerEnergy( int level) {

    

    return JAB_triplet( level  , level ) * votca::tools::conv::ryd2ev_f;
}



float BSECoupling::getTripletCouplingElement( int levelA, int levelB) {

    
   
    //cout << levelA<<endl;
    //cout << levelB + _levA<<endl;
    return JAB_triplet( levelA  , levelB + _levA ) * votca::tools::conv::ryd2ev_f;


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
bool BSECoupling::CalculateCouplings(Orbitals* _orbitalsA, Orbitals* _orbitalsB, Orbitals* _orbitalsAB) {
          
    LOG(logDEBUG,*_pLog) << "  Calculating exciton couplings" << flush;
    

    //check to see if ordering of atoms agrees
    const std::vector<QMAtom*> atomsA=_orbitalsA->QMAtoms();
    const std::vector<QMAtom*> atomsB=_orbitalsB->QMAtoms();
    const std::vector<QMAtom*> atomsAB=_orbitalsAB->QMAtoms();
    
    for (unsigned i=0;i<atomsAB.size();i++){
        QMAtom* dimer=atomsAB[i];
        QMAtom* monomer=NULL;
        if (i<atomsA.size()){
            monomer=atomsA[i];
        }
        else if (i<atomsB.size()+atomsA.size() ){
            monomer=atomsB[i-atomsA.size()];
        }
        else{
            throw runtime_error((boost::format("Number of Atoms in dimer %3i and the two monomers A:%3i B:%3i does not agree") %atomsAB.size() %atomsA.size() %atomsB.size()).str());
        }
        
        if(monomer->type != dimer->type){
            throw runtime_error("\nERROR: Atom types do not agree in dimer and monomers\n");
        }
        if(std::abs(monomer->x-dimer->x)>0.001 || std::abs(monomer->y-dimer->y)>0.001 || std::abs(monomer->z-dimer->z)>0.001){
            LOG(logERROR,*_pLog) << "======WARNING=======\n Coordinates of monomers and dimer atoms do not agree, do you know what you are doing?\n " << flush;
            break;
        }
        
    }
    
    
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
    
    if(_levA>_bseA_singlet_exc){
        LOG(logDEBUG,*_pLog) << "Number of excitons you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levA=_bseA_singlet_exc;
    }
    if(_levB>_bseB_singlet_exc){
        LOG(logDEBUG,*_pLog) << "Number of excitons you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levB=_bseB_singlet_exc;
    }
    
    if(_unoccA>_bseA_ctotal){
        LOG(logDEBUG,*_pLog) << "Number of occupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccA=_bseA_ctotal;
    }
    if(_unoccB>_bseB_ctotal){
        LOG(logDEBUG,*_pLog) << "Number of occupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccB=_bseB_ctotal;
    }
    if(_occA>_bseA_vtotal){
        LOG(logDEBUG,*_pLog) << "Number of unoccupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _occA=_bseA_vtotal;
    }
    if(_occB>_bseB_vtotal){
        LOG(logDEBUG,*_pLog) << "Number of unoccupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _occB=_bseB_vtotal;
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
    
    

    
    // DFT levels of monomers can be reduced to those used in BSE
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
    cout<< "_psi_AxB_dimer"<<endl;
    for (int i=0;i<_psi_AxB_dimer_basis.size1();i++){
        double mag=0.0;
        for (int j=0;j<_psi_AxB_dimer_basis.size2();j++){
            mag+=_psi_AxB_dimer_basis(i,j)*_psi_AxB_dimer_basis(i,j);
            
    }
         if (mag<0.95){
            throw runtime_error("\nERROR: Projection of monomer orbitals on dimer is insufficient, maybe the orbital order is screwed up, otherwise increase dimer basis.\n");
        }
    }
    //exit(0);
    // _psi_AxB_dimer_basis = T in notes, dimension ( LA + LB, LD)
    
    
   // cout << "Size of _psi_AxB_dimer_basis " << _psi_AxB_dimer_basis.size1() << " : " <<  _psi_AxB_dimer_basis.size2() << flush; 
    
    
    //notation AB is CT states with A+B-, BA is the counterpart
    //Setting up CT-states:
    LOG(logDEBUG,*_pLog) << "Setting up CT-states" << flush; 
    //Number of A+B- states
    int noAB=_occA*_unoccB;
    //Number of A-B+ states
    int noBA=_unoccA*_occB;
    
    
    ub::matrix<int> comb_CTAB;
    comb_CTAB.resize(noAB,2);
    _cnt = 0;
    

    
    
    // iterate A over occupied, B over unoccupied
    int v_start=_bseA_vtotal-_occA;
    for ( int _v = v_start; _v < _bseA_vtotal; _v++){
        for ( int _c = 0; _c <_unoccB; _c++){            
            comb_CTAB(_cnt,0) =_v;
            comb_CTAB(_cnt,1) = _bseA_vtotal+_bseA_ctotal+_bseB_vtotal + _c;
           
            _cnt++;
        }
    }
    LOG(logDEBUG,*_pLog) <<noAB <<" CT states A+B- created" << flush;
    //cout << "comb_CTAB" << endl;
    //cout << comb_CTAB << endl;
    
    ub::matrix<int> comb_CTBA;
    comb_CTBA.resize(noBA,2);
    _cnt = 0;
    // iterate A over unoccupied, B over occupied
    v_start=_bseB_vtotal-_occB;
    for ( int _v = v_start; _v < _bseB_vtotal; _v++){
        for ( int _c = 0; _c <_unoccA; _c++){            
            comb_CTBA(_cnt,0) =_bseA_vtotal+_bseA_ctotal+_v;
            comb_CTBA(_cnt,1) = _bseA_vtotal+ _c;
            
            _cnt++;
        }
    }
    LOG(logDEBUG,*_pLog) <<noBA <<" CT states B+A- created" << flush;
    //cout << "comb_CTBA" << endl;
    //cout << comb_CTBA << endl;
    
    
    
    // these 4 matrixes, matrix(i,j) contains the j-th dimer MO component of the i-th excitation
    ub::matrix<float> ctAB;
    ctAB.resize(noAB,_bseAB_size);
    for ( int _i_CT = 0 ; _i_CT < noAB ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctAB(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTAB(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTAB(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
    ub::matrix<float>ctBA;
    ctBA.resize(noBA,_bseAB_size);
    for ( int _i_CT = 0 ; _i_CT < noBA ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctBA(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTBA(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTBA(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
      
    // some more convenient storage
    
    ub::matrix<float> _kap;
    _kap.resize(_bseA_size,_bseAB_size);
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) );
            
        }
    }

    ub::matrix<float> check;
    check.resize(_bseA_size,_bseA_size);
    check= ub::prod(_kap,ub::trans(_kap));
    cout << "check"<<endl;
    
    cout <<ub::project(check, ub::range (0, 10 ), ub::range ( 0, 10 ) ) <<endl;
    
    
    ub::matrix<float> _kbp;
    _kbp.resize(_bseB_size,_bseAB_size);
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) );
        }
    }
    
    // Same routines but also take <v|c'> <c|v'> projections into account 
    /*
    ub::matrix<float> _kap;
    _kap.resize(_bseA_size,_bseAB_size);
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) )+
              _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,1) )     ;
            
        }
    }

    

    
    
    ub::matrix<float> _kbp;
    _kbp.resize(_bseB_size,_bseAB_size);
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) )+
                    _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,1) );
        }
    }
    */ 
  

    LOG(logDEBUG,*_pLog)  << "   construct projection of product functions " << flush; 

   
    
        
    //     cout << "Size of _kap " << _kap.size1() << " : " <<  _kap.size2() << "\n" << flush; 
    //     cout << "Size of _kbp " << _kbp.size1() << " : " <<  _kbp.size2() << "\n" << flush; 
    
    // now the different spin types
    if ( _doSinglets){
         LOG(logDEBUG,*_pLog)  << "   Evaluating singlets"  << flush; 
      // get singlet BSE Hamiltonian from _orbitalsAB
      ub::matrix<float> _Hamiltonian_AB = _eh_d + 2.0 * _eh_x;
        const ub::matrix<float>& _bseA = ub::project( _orbitalsA->BSESingletCoefficients(),
                ub::range (0, _orbitalsA->BSESingletCoefficients().size1() ), ub::range ( 0, _levA )  );
        const ub::matrix<float>& _bseB = ub::project( _orbitalsB->BSESingletCoefficients(),
                ub::range (0, _orbitalsB->BSESingletCoefficients().size1() ), ub::range ( 0, _levB )  );
        
      //  cout << "Size of _Hamiltonian_AB " << _Hamiltonian_AB.size1() << " : " <<  _Hamiltonian_AB.size2() << flush;
// cout << "Size of _bseA " << _bseA.size1() << " : " <<  _bseA.size2() << flush; 
   //      cout << "Size of _bseB " << _bseB.size1() << " : " <<  _bseB.size2() << flush; 

     
        bool _singlets = ProjectExcitons( _kap, _kbp,ctAB,ctBA, _bseA, _bseB, _Hamiltonian_AB, JAB_singlet);
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
    
                
        
    if ( _doTriplets){
        LOG(logDEBUG,*_pLog)  << "   Evaluating triplets" << flush; 
        // get triplet BSE Hamiltonian from _orbitalsAB
        ub::matrix<float> _Hamiltonian_AB = _eh_d;
        
     
        
        
        //ub::matrix<float>& _bseA= _orbitalsA->BSETripletCoefficients();
        
        
        const ub::matrix<float>& _bseA = ub::project( _orbitalsA->BSETripletCoefficients(),
                ub::range (0, _orbitalsA->BSETripletCoefficients().size1() ), ub::range ( 0, _levA )  );
        const ub::matrix<float>& _bseB = ub::project( _orbitalsB->BSETripletCoefficients(),
                ub::range (0, _orbitalsB->BSETripletCoefficients().size1() ), ub::range ( 0, _levB )  );
        //ub::matrix<float>& _bseB = _orbitalsB->BSETripletCoefficients();
        
  
        
        
        bool _triplets = ProjectExcitons( _kap, _kbp,ctAB,ctBA, _bseA, _bseB, _Hamiltonian_AB, JAB_triplet);
        if ( _triplets ) {
            LOG(logDEBUG,*_pLog)  << "   calculated triplet couplings " << flush;
        }
    }
    
    
    LOG(logDEBUG,*_pLog) << "  Done with exciton couplings" << flush;
    return true;   
};


bool BSECoupling::ProjectExcitons(const ub::matrix<float>& _kap,const ub::matrix<float>& _kbp,
                                  const ub::matrix<float>& ctAB,const ub::matrix<float>& ctBA, 
                                  const ub::matrix<float>& _bseA, const ub::matrix<float>& _bseB, 
                                  const ub::matrix<float>& _H, ub::matrix<float>& _J){
    
    
    //cout << " Dimensions of _bseA " << _bseA.size1() << " : " << _bseA.size2() << endl;
    
     // get projection of monomer excitons on dimer product functions
     ub::matrix<float> _proj_excA = ub::prod( ub::trans( _bseA ), _kap);
     ub::matrix<float> _proj_excB = ub::prod( ub::trans( _bseB ), _kbp);
     //cout << "_bseA"<<_bseA.size1()<<"x"<<_bseA.size2()<<endl;
     //cout << "_kap"<<_kap.size1()<<"x"<<_kap.size2()<<endl;  
     //cout << "_proj_excA"<<_proj_excA.size1()<<"x"<<_proj_excA.size2()<<endl;
     //cout << "_bseB"<<_bseB.size1()<<"x"<<_bseB.size2()<<endl;
     //cout << "_kbp"<<_kbp.size1()<<"x"<<_kbp.size2()<<endl;
     //cout << "_proj_excB"<<_proj_excB.size1()<<"x"<<_proj_excB.size2()<<endl;
     
     unsigned _bseA_exc = _bseA.size2();
     unsigned _bseB_exc = _bseB.size2();
     unsigned _bse_exc=_bseA_exc+_bseB_exc;
     unsigned _ctAB=ctAB.size1();
     
     unsigned _ctBA=ctBA.size1();
     unsigned _ct=_ctAB+_ctBA;
     //cout << _ctAB <<_ctBA<<endl;
     ub::matrix<float> _J_dimer = ub::zero_matrix<float>( _bse_exc +_ct, _bse_exc+_ct );
     ub::matrix<float> _S_dimer = ub::zero_matrix<float>( _bse_exc +_ct, _bse_exc +_ct);
     
     //cout << _J_dimer.size1()<< " : "<<_J_dimer.size2()<<endl;
     //cout << _S_dimer.size1()<< " : "<<_S_dimer.size2()<<endl;

      LOG(logDEBUG,*_pLog) << "Setting up coupling matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;
     // matrix _J
     
    //  E_A         J_AB        J_A_ABCT        J_A_BACT
    //  J_BA        E_B         J_B_ABCT        J_B_BACT
    //  J_ABCT_A    J_ABCT_B    E_ABCT          J_ABCT_BACT
    //  J_BACT_A   J_BACT_B    J_BACT_ABCT     E_BACT
     
     // I think this only works for hermitian/symmetric H so only in TDA
     // setup J
     
     // do not reorder ifs because of the temps, I would like to paralleize this part but i have no clue how to do it simply and elegantly
     
     ub::matrix<float> _temp = ub::prod( _H , ub::trans(_proj_excA) );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excA, _temp ); // E_A = proj_excA x H x trans(proj_excA)
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bse_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excB, _temp ); // J_BA = proj_excB x H x trans(proj_excA)
     if(_ctAB>0){
     ub::project( _J_dimer,  ub::range (_bse_exc, _bse_exc+_ctAB ), ub::range ( 0, _bseA_exc )  ) = ub::prod( ctAB, _temp );// J_ABCT_A
     }
     if(_ctBA>0){
     ub::project( _J_dimer,  ub::range ( _bse_exc+_ctAB, _bse_exc+_ct ), ub::range ( 0, _bseA_exc )  ) = ub::prod( ctBA, _temp );// J_BACT_A
     }
     
     
     _temp = ub::prod( _H , ub::trans(_proj_excB) );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bseA_exc, _bse_exc)  ) =ub::prod( _proj_excA, _temp ); // J_AB = proj_excA x H x trans(proj_excB)
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bse_exc ), ub::range ( _bseA_exc,_bse_exc)  ) = ub::prod(_proj_excB, _temp ); // E_B = proj_excB x H x trans(proj_excB)
     if(_ctAB>0){
     ub::project( _J_dimer,  ub::range (_bse_exc, _bse_exc+_ctAB ), ub::range ( _bseA_exc,_bse_exc)  ) = ub::prod(ctAB, _temp ); // J_ABCT_B
     }
     if(_ctBA>0){
     ub::project( _J_dimer,  ub::range (_bse_exc+_ctAB, _bse_exc+_ct ), ub::range ( _bseA_exc,_bse_exc)  ) = ub::prod(ctBA, _temp );// J_BACT_B
     }
     
     
     if(_ctAB>0){
      _temp = ub::prod( _H , ub::trans(ctAB) );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bse_exc, _bse_exc+_ctAB )  ) = ub::prod( _proj_excA, _temp );
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bse_exc ), ub::range ( _bse_exc, _bse_exc+_ctAB )  ) = ub::prod( _proj_excB, _temp );
     ub::project( _J_dimer,  ub::range (_bse_exc, _bse_exc+_ctAB ), ub::range ( _bse_exc, _bse_exc+_ctAB )  ) = ub::prod( ctAB, _temp );
     ub::project( _J_dimer,  ub::range (_bse_exc+_ctAB, _bse_exc+_ct ), ub::range ( _bse_exc, _bse_exc+_ctAB )  ) = ub::prod( ctBA, _temp );
     }
     
     
     if(_ctBA>0){
      _temp = ub::prod( _H , ub::trans(ctBA) );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bse_exc+_ctAB, _bse_exc+_ct )  ) = ub::prod( _proj_excA, _temp );
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bse_exc ), ub::range ( _bse_exc+_ctAB, _bse_exc+_ct )  ) = ub::prod( _proj_excB, _temp );
     ub::project( _J_dimer,  ub::range (_bse_exc, _bse_exc+_ctAB ), ub::range ( _bse_exc+_ctAB, _bse_exc+_ct )   ) = ub::prod( ctAB, _temp );
     ub::project( _J_dimer,  ub::range (_bse_exc+_ctAB, _bse_exc+_ct ), ub::range ( _bse_exc+_ctAB, _bse_exc+_ct )   ) = ub::prod( ctBA, _temp ); 
     }

     

    
    LOG(logDEBUG,*_pLog) << "Setting up overlap matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;
     // setup S
    
     _temp =ub::trans(_proj_excA);
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excA, _temp ); // E_A = proj_excA x H x trans(proj_excA)
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bse_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excB, _temp ); // J_BA = proj_excB x H x trans(proj_excA)
     if(_ctAB>0){
     ub::project( _S_dimer,  ub::range (_bse_exc, _bse_exc+_ctAB ), ub::range ( 0, _bseA_exc )  ) = ub::prod( ctAB, _temp );// J_ABCT_A
     }
     if(_ctBA>0){
     ub::project( _S_dimer,  ub::range ( _bse_exc+_ctAB, _bse_exc+_ct ), ub::range ( 0, _bseA_exc )  ) = ub::prod( ctBA, _temp );// J_BACT_A
     }
     
     
     _temp =  ub::trans(_proj_excB);
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bseA_exc, _bse_exc)  ) =ub::prod( _proj_excA, _temp ); // J_AB = proj_excA x H x trans(proj_excB)
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bse_exc ), ub::range ( _bseA_exc,_bse_exc)  ) = ub::prod(_proj_excB, _temp ); // E_B = proj_excB x H x trans(proj_excB)
     if(_ctAB>0){
     ub::project( _S_dimer,  ub::range (_bse_exc, _bse_exc+_ctAB ), ub::range ( _bseA_exc,_bse_exc)  ) = ub::prod(ctAB, _temp ); // J_ABCT_B
     }
     if(_ctBA>0){
     ub::project( _S_dimer,  ub::range (_bse_exc+_ctAB, _bse_exc+_ct ), ub::range ( _bseA_exc,_bse_exc)  ) = ub::prod(ctBA, _temp );// J_BACT_B
     }
     
     if(_ctAB>0){
      _temp =  ub::trans(ctAB);
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bse_exc, _bse_exc+_ctAB )  ) = ub::prod( _proj_excA, _temp );
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bse_exc ), ub::range ( _bse_exc, _bse_exc+_ctAB )  ) = ub::prod( _proj_excB, _temp );
     ub::project( _S_dimer,  ub::range (_bse_exc, _bse_exc+_ctAB ), ub::range ( _bse_exc, _bse_exc+_ctAB )  ) = ub::prod( ctAB, _temp );
     ub::project( _S_dimer,  ub::range (_bse_exc+_ctAB, _bse_exc+_ct ), ub::range ( _bse_exc, _bse_exc+_ctAB )  ) = ub::prod( ctBA, _temp );
     }
     if(_ctBA>0){
      _temp = ub::trans(ctBA);
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bse_exc+_ctAB, _bse_exc+_ct )  ) = ub::prod( _proj_excA, _temp );
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bse_exc ), ub::range ( _bse_exc+_ctAB, _bse_exc+_ct )  ) = ub::prod( _proj_excB, _temp );
     ub::project( _S_dimer,  ub::range (_bse_exc, _bse_exc+_ctAB ), ub::range ( _bse_exc+_ctAB, _bse_exc+_ct )   ) = ub::prod( ctAB, _temp );
     ub::project( _S_dimer,  ub::range (_bse_exc+_ctAB, _bse_exc+_ct ), ub::range ( _bse_exc+_ctAB, _bse_exc+_ct )   ) = ub::prod( ctBA, _temp ); 
     }
    
     // cout<<_S_dimer<<endl;
     ub::vector<float> _S_eigenvalues; 
     cout <<endl;
      cout<<"_J_dimer"<<endl;
     cout << _J_dimer*conv::ryd2ev_f<<endl;
     cout << "S"<<endl;
     cout << _S_dimer<< endl;
     linalg_eigenvalues( _S_eigenvalues, _S_dimer);
      
     LOG(logDEBUG,*_pLog) << "Smallest value of dimer overlapmatrix is "<< _S_eigenvalues[0]<< flush;

     if ( _S_eigenvalues[0] < 0.0 ) {
         
         cerr << " \n Negative eigenvalues in projected overlap. Projection quality not sufficient. Try increasing dimer basis. " << endl;
         return false;
         
     }
     
     ub::matrix<float> _diagS = ub::zero_matrix<float>( _bse_exc + _ct  , _bse_exc + _ct );
     for ( unsigned _i =0; _i < _bse_exc + _ct ; _i++){

         _diagS(_i,_i) = 1.0/sqrt(_S_eigenvalues[_i]);
     }
     
     //ub::matrix<float> _transform = ub::prod( _S_dimer, ub::prod( _diagS, ub::trans(_S_dimer) )  );
     //cout << endl;
     //cout << "J_dimer"<<endl;
     //cout <<_J_dimer<<endl;
     ub::matrix<float> _transtemp = ub::prod( _diagS, ub::trans(_S_dimer));
     ub::matrix<float> _transform = ub::prod( _S_dimer,_transtemp );

     ub::matrix<float> _J_temp = ub::prod(_J_dimer, _transform);
     
             // final coupling elements
     // _J = ub::prod( _transform, ub::prod(_J_dimer, _transform));
    ub::matrix<float>_J_ortho = ub::prod( _transform, _J_temp);
    cout<<endl;
    cout<< "diagonals"<<endl;
    for (int i=0;i<_J_ortho.size1();i++){
        cout << _J_ortho(i,i)*conv::ryd2ev_f<<" ";
    }
    cout <<endl;
    cout << "J_ortho"<<endl;
    cout << _J_ortho*conv::ryd2ev_f<<endl;
     //Setting up effective Coupling matrix
     LOG(logDEBUG,*_pLog) << "Setting up effective/perturbative Coupling matrix of size: "<< _bse_exc  <<"x"<<  _bse_exc << flush;
     
     
     
     
     ub::vector<float> _J_eigenvalues;
     linalg_eigenvalues(_J_eigenvalues,_J_ortho);
     
     //finding the eigenstates which are closest to the the original states
     
     cout << "_J_eigenvalues" << endl;
    cout << _J_eigenvalues*conv::ryd2ev_f << endl;
     cout << "_J_eigenvectors" << endl;
     cout << _J_ortho<<endl;
     
     
     //setting up transformation matrix _T and diagonal matrix _E for the eigenvalues;
     
     ub::matrix<float> _E=ub::zero_matrix<float>(_bse_exc,_bse_exc);
     ub::matrix<float> _T=ub::zero_matrix<float>(_bse_exc,_bse_exc);
     //find the eigenvectors which are most similar to the initial states
     
     LOG(logDEBUG,*_pLog) << "Sorting states according to similiarity with the FE states " << flush;
     
     std::vector<unsigned> index;
     //column
      for (unsigned i = 0; i < _bse_exc; i++) {
                float close = 0.0;
                unsigned ind = 0;
                //row
                for (unsigned j = 0; j < _bse_exc + _ct; j++) {
                    bool check=true;
                    // if index i is already in index
                    // should not happen but if one vector was similar to tow others.
                    for (unsigned l=0;l<index.size();l++){
                        if (j==index[l]){
                            check=false;
                            break;
                        }
                    }
                    
                    if (check && std::abs(_J_ortho(i, j)) > close) {
                        ind = j;
                        close=std::abs(_J_ortho(i, j));
                    }
                }
                index.push_back(ind);
            }
     
     LOG(logDEBUG,*_pLog) << "Order is: [Initial state n->nth eigenvalue]"<<flush;
     for (unsigned i=0;i<index.size();i++){
         if(i<_bseA_exc){
      LOG(logDEBUG,*_pLog) <<"    A"<<i+1<<":"<<i+1<<"->"<<index[i]+1<<" " ;   
         }
         else{
     LOG(logDEBUG,*_pLog) <<"    B"<<i+1-_bseA_exc<<":"<<i+1<<"->"<<index[i]+1<<" " ;  
         }
                 
     }
     LOG(logDEBUG,*_pLog)<< flush;
     //row
     for (unsigned i = 0; i < _bse_exc; i++) {
         unsigned k=index[i];
                float norm = 0.0;
                //column
                for (unsigned j = 0; j < _bse_exc; j++) {
                    
                    norm += _J_ortho(j, k)*_J_ortho(j, k);
                }
                for (unsigned j = 0; j < _bse_exc; j++) {
                    if (i == j) {
                        _E(i, i) = _J_eigenvalues(k);
                    }
                    _T(j,i ) = _J_ortho(j,k) / std::sqrt(norm);
                }
            }
     cout << "_E" <<endl;
     cout << _E*conv::ryd2ev_f <<endl;
          cout << "_T" <<endl;

     cout << _T << endl;
     _temp=ub::prod(_E,ub::trans(_T));
     cout << "_J" <<endl;
     _J=ub::prod(_T,_temp);
     cout <<_J*conv::ryd2ev_f<<endl;
     
     /*
     for (unsigned i=0;i<_J.size1();i++){
         for (unsigned j=0;j<i;j++){
             float epsilonA=_J_ortho(i,i);
             float epsilonB=_J_ortho(j,j);
           _J(i,j)=_J_ortho(i,j);
           for(int k=_bse_exc;k<_bse_exc+_ct;k++){
               float omegaAB=_J_ortho(k,k);
             _J(i,j)+= _J_ortho(i,k)* _J_ortho(j,k)*(1/(epsilonA-omegaAB)+1/(epsilonB-omegaAB));
           }
          _J(j,i)=_J(i,j);
                 
                 
           
        }
     }
      for (unsigned i=0;i<_J.size1();i++){
          _J(i,i)=_J_ortho(i,i);
      }
     */ 
     
    
   
     /*
     // This is not necessary right now, as we use perturbation theory on the non-orthogonal system to derive the couplings 
     // see [1]B. Lunkenheimer, “Simulationen zur Exzitonendiffusion in organischen Halbleitern,” Universitätsbibliothek Mainz, 2014.

     ub::matrix<float> _diagS = ub::zero_matrix<float>( _bse_exc + _ct  , _bse_exc + _ct );
     for ( int _i =0; _i < _bse_exc + _ct ; _i++){

         _diagS(_i,_i) = 1.0/sqrt(_S_eigenvalues[_i]);
     }
     
     //ub::matrix<float> _transform = ub::prod( _S_dimer, ub::prod( _diagS, ub::trans(_S_dimer) )  );

     ub::matrix<float> _transtemp = ub::prod( _diagS, ub::trans(_S_dimer));
     ub::matrix<float> _transform = ub::prod( _S_dimer,_transtemp );

     ub::matrix<float> _J_temp = ub::prod(_J_dimer, _transform);
     
             // final coupling elements
     // _J = ub::prod( _transform, ub::prod(_J_dimer, _transform));
    _J = ub::prod( _transform, _J_temp);
     */
      
     return true;
    
    
    
    
}
    
}}
