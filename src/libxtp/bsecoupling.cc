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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>

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
namespace CTP = votca::ctp;
using boost::format;

void BSECoupling::Initialize(Property* options){
    
    #if (GWBSE_DOUBLE)
        LOG(CTP::logDEBUG, *_pLog) <<  " Compiled with full double support" << flush;   
    #else
        LOG(CTP::logDEBUG, *_pLog) <<  " Compiled with float/double mixture (standard)" << flush;   
    #endif
    
    std::string key = Identify(); 
    _doSinglets=false;
    _doTriplets=false;
    _do_perturbation=false;
    _do_full_diag=false;
    
    _openmp_threads = 0;
    
    if ( options->exists( key + ".openmp") ) {
                 _openmp_threads = options->get(key + ".openmp").as<int> ();
            }
    
    string spintype   = options->get(key + ".spin").as<string> ();
        if(spintype=="all"){
            _doSinglets=true;
            _doTriplets=true;
        }
        else if(spintype=="triplet"){
            _doTriplets=true;
        }
        else if(spintype=="singlet"){
            _doSinglets=true;
        }
        else{
            throw std::runtime_error((boost::format("Choice % for type not known. Available singlet,triplet,all") % spintype).str());
        }
    _degeneracy = options->get(key + ".degeneracy").as<double> ();
     string algorithm;
   if ( options->exists( key + ".algorithm") ) {
                algorithm = options->get(key + ".algorithm").as<string> ();
                 if(algorithm=="perturbation"){
                    _do_perturbation=true;
           
                }
                else if(algorithm=="fulldiag"){
                    _do_full_diag=true;
            }
        else{
            throw std::runtime_error((boost::format("Choice % for algorithm not known. Available fulldiag(default),perturbation.") % algorithm).str());
        }
            }
   else{
       _do_full_diag=true;
       algorithm="fulldiag";
   }
    
    LOG(CTP::logDEBUG, *_pLog) <<  " Using "<< algorithm << " to solve for couplings"<< flush;   
    
        _levA  = options->get(key + ".moleculeA.states").as<int> ();
        _levB  = options->get(key + ".moleculeB.states").as<int> ();
        _occA  = options->get(key + ".moleculeA.occLevels").as<int> ();
        _occB  = options->get(key + ".moleculeB.occLevels").as<int> ();
        _unoccA  = options->get(key + ".moleculeA.unoccLevels").as<int> ();
        _unoccB  = options->get(key + ".moleculeB.unoccLevels").as<int> ();
        
        if ( options->exists(key+".moleculeA.Frenkel")) {
         _FeA = options->get(key + ".moleculeA.Frenkel").as<int> ();
        
        }
        else{
            _FeA=_levA;
        }
        if ( options->exists(key+".moleculeB.Frenkel")) {
         _FeB = options->get(key + ".moleculeB.Frenkel").as<int> ();
        
        }
        else{
            _FeB=_levB;
        }
        
        
}

void BSECoupling::addoutput(Property *_type_summary,Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB){
    //cout << JAB_singlet<<endl;
    //cout << JAB_singlet*conv::ryd2ev_f<<endl;
    string algorithm="";
    if (_do_perturbation){
        algorithm="perturbation";
    }
    else if(_do_full_diag){
        algorithm="full_diag";
    }
    _type_summary->setAttribute("algorithm",algorithm);
    if (_doSinglets){
        Property *_singlet_summary = &_type_summary->add("singlets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB <_levB ; ++stateB ) {
               real_gwbse JAB = getSingletCouplingElement( stateA , stateB );
               //real_gwbse_gwbse energyAD = getSingletDimerEnergy( stateA  );
               //real_gwbse energyBD = getSingletDimerEnergy( stateB  );
               Property *_coupling_summary = &_singlet_summary->add("coupling", boost::lexical_cast<string>(JAB)); 
               real_gwbse energyA = _orbitalsA->BSESingletEnergies()(stateA)*conv::ryd2ev_f;
               real_gwbse energyB = _orbitalsB->BSESingletEnergies()(stateB)*conv::ryd2ev_f;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", energyA);
               _coupling_summary->setAttribute("energyB", energyB);
               //_coupling_summary->setAttribute("energyA_dimer", energyAD);
               //_coupling_summary->setAttribute("energyB_dimer", energyBD);
           } 
        }
    }
    
    //cout << JAB_triplet<<endl;
    //cout << JAB_triplet*conv::ryd2ev_f<<endl;
    if ( _doTriplets){
        Property *_triplet_summary = &_type_summary->add("triplets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB < _levA ; ++stateB ) {
               real_gwbse JAB = getTripletCouplingElement( stateA , stateB );
               //real_gwbse energyAD = getTripletDimerEnergy( stateA  );
               //real_gwbse energyBD = getTripletDimerEnergy( stateB  );
               Property *_coupling_summary = &_triplet_summary->add("coupling", boost::lexical_cast<string>(JAB)); 
               real_gwbse energyA = _orbitalsA->BSETripletEnergies()(stateA)*conv::ryd2ev_f;
               real_gwbse energyB = _orbitalsB->BSETripletEnergies()(stateB)*conv::ryd2ev_f;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", energyA);
               _coupling_summary->setAttribute("energyB", energyB);
               //_coupling_summary->setAttribute("energyA_dimer", energyAD);
               //_coupling_summary->setAttribute("energyB_dimer", energyBD);
           } 
        }
    }       
}


real_gwbse BSECoupling::getSingletCouplingElement( int levelA, int levelB) {

    

    return JAB_singlet( levelA  , levelB +  _levA ) * votca::tools::conv::ryd2ev_f;
}
real_gwbse BSECoupling::getSingletDimerEnergy( int level) {

    

    return JAB_singlet( level  , level ) * votca::tools::conv::ryd2ev_f;
}
real_gwbse BSECoupling::getTripletDimerEnergy( int level) {

    

    return JAB_triplet( level  , level ) * votca::tools::conv::ryd2ev_f;
}



real_gwbse BSECoupling::getTripletCouplingElement( int levelA, int levelB) {

    
   
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
    Orbitals* _orbitalsAB, ub::matrix<real_gwbse>* _JAB) {

    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "Calculating exciton couplings" << flush;
    


        
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        LOG(CTP::logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
        return false;
    }

    // number of levels stored in monomers
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    boost::timer t; // start timing
    double _st = t.elapsed();
    
    
    
        
    // get exciton information of molecule A
    ub::matrix<real_gwbse>& _bseA = _orbitalsA->BSESingletCoefficients();
    int _bseA_cmax = _orbitalsA->getBSEcmax();
    int _bseA_cmin = _orbitalsA->getBSEcmin();
    int _bseA_vmax = _orbitalsA->getBSEvmax();
    int _bseA_vmin = _orbitalsA->getBSEvmin();
    int _bseA_vtotal = _bseA_vmax - _bseA_vmin +1 ;
    int _bseA_ctotal = _bseA_cmax - _bseA_cmin +1 ;
    int _bseA_size   = _bseA_vtotal * _bseA_ctotal;
    int _bseA_exc    = _bseA.size2();
    
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " molecule A has " << _bseA_exc << " excitons with dimension " << _bseA_size << flush;
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
    ub::matrix<real_gwbse>& _bseB = _orbitalsB->BSESingletCoefficients();
    int _bseB_cmax = _orbitalsB->getBSEcmax();
    int _bseB_cmin = _orbitalsB->getBSEcmin();
    int _bseB_vmax = _orbitalsB->getBSEvmax();
    int _bseB_vmin = _orbitalsB->getBSEvmin();
    int _bseB_vtotal = _bseB_vmax - _bseB_vmin +1 ;
    int _bseB_ctotal = _bseB_cmax - _bseB_cmin +1 ;
    int _bseB_size   = _bseB_vtotal * _bseB_ctotal;
    int _bseB_exc    = _bseB.size2();
    
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " molecule B has " << _bseB_exc << " excitons with dimension " << _bseB_size << flush;
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
    ub::matrix<real_gwbse>& _bseAB = _orbitalsAB->BSESingletCoefficients();
    int _bseAB_cmax = _orbitalsAB->getBSEcmax();
    int _bseAB_cmin = _orbitalsAB->getBSEcmin();
    int _bseAB_vmax = _orbitalsAB->getBSEvmax();
    int _bseAB_vmin = _orbitalsAB->getBSEvmin();
    int _bseAB_vtotal = _bseAB_vmax - _bseAB_vmin +1 ;
    int _bseAB_ctotal = _bseAB_cmax - _bseAB_cmin +1 ;
    int _bseAB_size   = _bseAB_vtotal * _bseAB_ctotal;
    int _bseAB_exc    = _bseAB.size2();
    
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " dimer AB has " << _bseAB_exc << " excitons with dimension " << _bseAB_size << flush;
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
    
    
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "Levels stored: Basis A[" << _levelsA << ":" << _basisA << "]"
                                     << " B[" << _levelsB << ":" << _basisB << "]" << flush;
    
    // DFT levels can be reduced to those used in BSE
    _levelsA = _bseA_vtotal + _bseA_ctotal;
    _levelsB = _bseB_vtotal + _bseB_ctotal;
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << " Levels used in BSE of molA: " << _bseA_vmin << " to " << _bseA_cmax << " total: " << _bseA_vtotal + _bseA_ctotal <<  flush;
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << " Levels used in BSE of molB: " << _bseB_vmin << " to " << _bseB_cmax << " total: " << _bseB_vtotal + _bseB_ctotal <<  flush;
    
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        LOG(CTP::logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
        return false;
    } 
    
    //       | Orbitals_A          0 |      | Overlap_A |     
    //       | 0          Orbitals_B |  X   | Overlap_B |  X  Transpose( Orbitals_AB )
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
    ub::matrix<double> _psi_AxB ( _levelsA + _levelsB, _basisA + _basisB  );
    
     
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "Constructing direct product AxB [" 
            << _psi_AxB.size1() << "x" 
            << _psi_AxB.size2() << "]";    
    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    
     
    //ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = *_orbitalsA->getOrbitals();
    // ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = *_orbitalsB->getOrbitals(); 

    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = ub::project( *_orbitalsA->getOrbitals() , ub::range(_bseA_vmin, _bseA_cmax+1) , ub::range ( 0, _basisA ));
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = ub::project( *_orbitalsB->getOrbitals(), ub::range(_bseB_vmin, _bseB_cmax+1) , ub::range ( 0, _basisB )); 
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();

    
    // psi_AxB * S_AB * psi_AB
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "Projecting  monomer orbitals onto dimer "; 
    ub::matrix<double> _orbitalsAB_Transposed = ub::trans( *_orbitalsAB->getOrbitals() );  
    if ( (*_orbitalsAB->getOverlap()).size1() == 0 ) {
            LOG(CTP::logERROR,*_pLog) << "Overlap matrix is not stored"; 
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
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();    

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
 
        
        //cout << " Monomer level projection (occ:virt:total)" << i << " : " << _proj_check(i,0) << " : " << _proj_check(i,1) << " : " << _proj_check(i,0) + _proj_check(i,1) << endl; 
        
        
    }
    
    
 
    
    
    
    
    
    
    
    
    
    // some more convenient storage
    
    ub::matrix<real_gwbse> _kap;
    _kap.resize(_bseA_size,_bseAB_size);
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) );
            
        }
    }

    
    ub::matrix<real_gwbse> _kbp;
    _kbp.resize(_bseB_size,_bseAB_size);
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) );
        }
    }

    
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " projection of product functions(" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();    
    
     
     // now project the exciton functions (ugly!)
     // ub::matrix<real_gwbse> _proj_excA;
     // _proj_excA.resize(_bseA_exc,_bseAB_exc);
     ub::matrix<real_gwbse> _temp = ub::prod( _kap, _bseAB );
     ub::matrix<real_gwbse> _proj_excA = ub::prod( ub::trans(_bseA), _temp);
     
     /*
     for ( unsigned i = 0 ; i < _proj_excA.size1() ; i++ ){
     
        real_gwbse check = 0.0;
        
        for ( unsigned j = 0; j < _proj_excA.size2(); j++ ){
            
            check += _proj_excA(i,j) * _proj_excA(i,j) ;
        }
        
        cout << " Monomer exciton " << i << " norm " << check << endl;
        
     }
     
     */
     
     _temp = ub::prod( _kbp, _bseAB );
     ub::matrix<real_gwbse> _proj_excB = ub::prod( ub::trans(_bseB ), _temp);
     
     
    // cout << "exciton projectors: " << _proj_excA.size1() << " : " << _proj_excA.size2();
     
     
     
     
     
     
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " projection of exciton wave functions (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();  
     
     
     ub::matrix<real_gwbse> _Hamiltonian_AB = ub::zero_matrix<real_gwbse>(_bseAB_size,_bseAB_size );
     ub::vector<real_gwbse>& _bseAB_energies = _orbitalsAB->BSESingletEnergies();

     for ( int _i=0; _i<_bseAB_size; _i++){
         _Hamiltonian_AB(_i,_i) = _bseAB_energies(_i);
     }
     
 
     
     ub::matrix<real_gwbse> _J_dimer = ub::zero_matrix<real_gwbse>( _bseA_exc + _bseB_exc , _bseA_exc + _bseB_exc );
     ub::matrix<real_gwbse> _S_dimer = ub::zero_matrix<real_gwbse>( _bseA_exc + _bseB_exc , _bseA_exc + _bseB_exc );
     
     // setup J
     _temp = ub::prod( _Hamiltonian_AB , ub::trans(_proj_excA) );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excA, _temp ); // J_AA = proj_excA x H x trans(proj_excA)
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excB, _temp ); // J_BA = proj_excB x H x trans(proj_excA)
     
     _temp = ub::prod( _Hamiltonian_AB , ub::trans(_proj_excB) );
     ub::project( _J_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc )  ) = ub::prod( _proj_excA, _temp ); // J_AB = proj_excA x H x trans(proj_excB)
     ub::project( _J_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc)  ) = ub::prod( _proj_excB, _temp ); // J_BB = proj_excB x H x trans(proj_excB)
     
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " setup J (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();  

     
     // setup S
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excA, ub::trans(_proj_excA) ); // S_AA = proj_excA x trans(proj_excA)
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( 0, _bseA_exc )  ) = ub::prod( _proj_excB, ub::trans(_proj_excA)); // S_BA = proj_excB x trans(proj_excA)
     ub::project( _S_dimer,  ub::range (0, _bseA_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc )  ) = ub::prod( _proj_excA, ub::trans(_proj_excB)); // J_AB = proj_excA x H x trans(proj_excB)
     ub::project( _S_dimer,  ub::range (_bseA_exc, _bseA_exc + _bseB_exc ), ub::range ( _bseA_exc, _bseA_exc + _bseB_exc)  ) = ub::prod( _proj_excB, ub::trans(_proj_excB) ); // J_BB = proj_excB x H x trans(proj_excB)
     
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " setup S (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();  
     
     ub::vector<real_gwbse> _S_eigenvalues; 
     linalg_eigenvalues( _S_eigenvalues, _S_dimer);

     ub::matrix<real_gwbse> _diagS = ub::zero_matrix<real_gwbse>( _bseA_exc + _bseB_exc  , _bseA_exc + _bseB_exc );
     for ( int _i =0; _i < _bseA_exc + _bseB_exc ; _i++){
         _diagS(_i,_i) = 1.0/sqrt(_S_eigenvalues[_i]);
     }
     _temp=ub::prod( _diagS, ub::trans(_S_dimer) ) ;
     ub::matrix<real_gwbse> _transform = ub::prod( _S_dimer, _temp );
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " transformation matrix (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();  
     // final coupling elements
     _temp=ub::prod(_J_dimer, _transform);
     ub::matrix<real_gwbse> _J_eff = ub::prod( _transform, _temp);
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " effective couplings (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed();
     
     
     //LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " singlet coupling: " << _J_eff(0,_bseA_exc+1)*13.6058 << " and " <<  _J_eff(_bseA_exc+1, 0) * 13.6058 << endl;
     //LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " singlet coupling: " << _J_eff(0,_bseA_exc)*13.6058 << " and " <<  _J_eff(_bseA_exc, 0) * 13.6058 << endl; 
     //LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " singlet coupling: " << _J_eff(1,_bseA_exc+1)*13.6058 << " and " <<  _J_eff(_bseA_exc+1, 1) * 13.6058 << endl; 


  
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " (" << t.elapsed() - _st << "s) " << flush; _st = t.elapsed(); 
    
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "Done with exciton couplings" << flush;
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
       LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Calculating exciton couplings" << flush;
     // set the parallelization 
    #ifdef _OPENMP
    
    if ( _openmp_threads > 0 ) omp_set_num_threads(_openmp_threads);      
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << " Using "<< omp_get_max_threads()<<" threads" << flush;
    #endif
    
    
    
    _orbitalsAB->setCoupledExcitonsA(_levA);
    _orbitalsAB->setCoupledExcitonsB(_levB);
    //check to see if ordering of atoms agrees
    const std::vector<CTP::QMAtom*> atomsA=_orbitalsA->QMAtoms();
    const std::vector<CTP::QMAtom*> atomsB=_orbitalsB->QMAtoms();
    const std::vector<CTP::QMAtom*> atomsAB=_orbitalsAB->QMAtoms();
    
    for (unsigned i=0;i<atomsAB.size();i++){
        CTP::QMAtom* dimer=atomsAB[i];
        CTP::QMAtom* monomer=NULL;
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
            LOG(CTP::logERROR,*_pLog) << "======WARNING=======\n Coordinates of monomers and dimer atoms do not agree, do you know what you are doing?\n " << flush;
            break;
        }
        
    }
    
    
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        LOG(CTP::logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
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

    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   molecule A has " << _bseA_singlet_exc << " singlet excitons with dimension " << _bseA_size << flush;
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   molecule A has " << _bseA_triplet_exc << " triplet excitons with dimension " << _bseA_size << flush;
    
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
    

    
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   molecule B has " << _bseB_singlet_exc << " singlet excitons with dimension " << _bseB_size << flush;
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   molecule B has " << _bseB_triplet_exc << " triplet excitons with dimension " << _bseB_size << flush;
    
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
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Number of excitons you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levA=_bseA_singlet_exc;
    }
    if(_levB>_bseB_singlet_exc){
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Number of excitons you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levB=_bseB_singlet_exc;
    }
    
    
    if(_FeA>_bseA_singlet_exc){
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Number of Frenkel states you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _FeA=_bseA_singlet_exc;
    }
    if(_FeB>_bseB_singlet_exc){
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Number of Frenkel states you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _FeB=_bseB_singlet_exc;
    }
    
    if(_unoccA>_bseA_ctotal){
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Number of occupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccA=_bseA_ctotal;
    }
    else if (_unoccA<0){
        _unoccA=_bseA_ctotal;
    }
    if(_unoccB>_bseB_ctotal){
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Number of occupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccB=_bseB_ctotal;
    }
    else if (_unoccB<0){
        _unoccB=_bseB_ctotal;
    }
    
    if(_occA>_bseA_vtotal){
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Number of unoccupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _occA=_bseA_vtotal;
    }
    else if (_occA<0){
        _occA=_bseA_vtotal;
    }
    if(_occB>_bseB_vtotal){
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Number of unoccupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _occB=_bseB_vtotal;
    }else if (_occB<0){
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
        LOG(CTP::logERROR,*_pLog) << "BSE EH int not stored in dimer " << flush;
        return false;
    }
    ub::matrix<real_gwbse>&    _eh_d = _orbitalsAB->eh_d(); 
    ub::matrix<real_gwbse>&    _eh_x = _orbitalsAB->eh_x(); 
    
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   dimer AB has BSE EH interaction (direct)   with dimension " << _eh_d.size1() << " x " <<  _eh_d.size2() << flush;
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   dimer AB has BSE EH interaction (exchange) with dimension " << _eh_x.size1() << " x " <<  _eh_x.size2() << flush;
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
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   levels used in BSE of molA: " << _bseA_vmin << " to " << _bseA_cmax << " total: " << _bseA_vtotal + _bseA_ctotal <<  flush;
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   levels used in BSE of molB: " << _bseB_vmin << " to " << _bseB_cmax << " total: " << _bseB_vtotal + _bseB_ctotal <<  flush;
    
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        LOG(CTP::logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
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
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   projecting monomer onto dimer orbitals" << flush; 
    ub::matrix<double> _orbitalsAB_Transposed = ub::trans( *_orbitalsAB->getOrbitals() );  
    if ( (*_orbitalsAB->getOverlap()).size1() == 0 ) {
            LOG(CTP::logERROR,*_pLog) << "Overlap matrix is not stored"; 
            return false;
    }
   
    ub::matrix<double> _psi_AB = ub::prod( *_orbitalsAB->getOverlap(), _orbitalsAB_Transposed );  
    ub::matrix<double> _psi_AxB_dimer_basis = ub::prod( _psi_AxB, _psi_AB );  
    _psi_AB.clear();
    //cout<< "_psi_AxB_dimer"<<endl;
    unsigned int LevelsA = _levelsA;
    for (unsigned i=0;i<_psi_AxB_dimer_basis.size1();i++){
        double mag=0.0;
        for (unsigned j=0;j<_psi_AxB_dimer_basis.size2();j++){
            mag+=_psi_AxB_dimer_basis(i,j)*_psi_AxB_dimer_basis(i,j);
            
    }
         if (mag<0.95){
            int monomer = 0;
            int level = 0;
            if ( i < LevelsA ) {
                monomer = 1;
                level   = _bseA_vmin + i;
            } else {
                monomer = 2;
                level   = _bseB_vmin + i -_levelsA;
                
            }
            LOG(CTP::logERROR,*_pLog) << "\nERROR: " << i << " Projection of orbital " << level << " of monomer " << monomer << " on dimer is insufficient,mag="<<mag<<" maybe the orbital order is screwed up, otherwise increase dimer basis.\n"<<flush;
        }
    }
    //exit(0);
    // _psi_AxB_dimer_basis = T in notes, dimension ( LA + LB, LD)
    
    
   // cout << "Size of _psi_AxB_dimer_basis " << _psi_AxB_dimer_basis.size1() << " : " <<  _psi_AxB_dimer_basis.size2() << flush; 
    
    
    //notation AB is CT states with A+B-, BA is the counterpart
    //Setting up CT-states:
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   Setting up CT-states" << flush; 
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
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  <<"  "<<noAB <<" CT states A+B- created" << flush;
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
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  <<"  "<<noBA <<" CT states B+A- created" << flush;
    //cout << "comb_CTBA" << endl;
    //cout << comb_CTBA << endl;
    
    
    
    // these 4 matrixes, matrix(i,j) contains the j-th dimer MO component of the i-th excitation
    ub::matrix<real_gwbse> ctAB;
    ctAB.resize(noAB,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_CT = 0 ; _i_CT < noAB ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctAB(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTAB(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTAB(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
    ub::matrix<real_gwbse>ctBA;
    ctBA.resize(noBA,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_CT = 0 ; _i_CT < noBA ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctBA(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTBA(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTBA(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
      
    // some more convenient storage
    
    ub::matrix<real_gwbse> _kap;
    _kap.resize(_bseA_size,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) );
            
        }
    }

    //ub::matrix<real_gwbse> check;
    //check.resize(_bseA_size,_bseA_size);
    //check= ub::prod(_kap,ub::trans(_kap));
    //cout << "check"<<endl;
    
    //cout <<ub::project(check, ub::range (0, 10 ), ub::range ( 0, 10 ) ) <<endl;
    
    
    ub::matrix<real_gwbse> _kbp;
    _kbp.resize(_bseB_size,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) );
        }
    }
    
    // Same routines but also take <v|c'> <c|v'> projections into account 
    /*
    ub::matrix<real_gwbse> _kap;
    _kap.resize(_bseA_size,_bseAB_size);
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) )+
              _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,1) )     ;
            
        }
    }

    

    
    
    ub::matrix<real_gwbse> _kbp;
    _kbp.resize(_bseB_size,_bseAB_size);
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) )+
                    _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,1) );
        }
    }
    */ 
  

    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   construct projection of product functions " << flush; 

   
    
        
    //     cout << "Size of _kap " << _kap.size1() << " : " <<  _kap.size2() << "\n" << flush; 
    //     cout << "Size of _kbp " << _kbp.size1() << " : " <<  _kbp.size2() << "\n" << flush; 
    
    // now the different spin types
    if ( _doSinglets){
         LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   Evaluating singlets"  << flush; 
      // get singlet BSE Hamiltonian from _orbitalsAB
      ub::matrix<real_gwbse> _Hamiltonian_AB = _eh_d + 2.0 * _eh_x;
        const ub::matrix<real_gwbse>& _bseA = ub::project( _orbitalsA->BSESingletCoefficients(),
                ub::range (0, _orbitalsA->BSESingletCoefficients().size1() ), ub::range ( 0, _FeA )  );
        const ub::matrix<real_gwbse>& _bseB = ub::project( _orbitalsB->BSESingletCoefficients(),
                ub::range (0, _orbitalsB->BSESingletCoefficients().size1() ), ub::range ( 0, _FeB )  );
        
      //  cout << "Size of _Hamiltonian_AB " << _Hamiltonian_AB.size1() << " : " <<  _Hamiltonian_AB.size2() << flush;
// cout << "Size of _bseA " << _bseA.size1() << " : " <<  _bseA.size2() << flush; 
   //      cout << "Size of _bseB " << _bseB.size1() << " : " <<  _bseB.size2() << flush; 

     
        bool _singlets = ProjectExcitons( _kap, _kbp,ctAB,ctBA, _bseA, _bseB, _Hamiltonian_AB, JAB_singlet);
        if ( _singlets ) {
            LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   calculated singlet couplings " << flush;
        }
        /*
        // diagonalize the effective Hamiltonian?
        ub::vector<real_gwbse> _coupled_energies;
        ub::matrix<real_gwbse> _coupled_coefficients;
	std::vector< std::vector<double> >& _free_dipolesA = _orbitalsA->TransitionDipoles();
	std::vector< std::vector<double> >& _free_dipolesB = _orbitalsB->TransitionDipoles();
        linalg_eigenvalues(*_JAB_singlet, _coupled_energies, _coupled_coefficients,_JAB_singlet->size1() );
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   calculated EVs of coupling matrix " << flush;
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
        //LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << " singlet coupling: " << _JAB_singlet->at_element(0,_bseA_singlet_exc)*13.6058 << " and " <<  _JAB_singlet->at_element(_bseA_singlet_exc, 0) * 13.6058 << endl; 
	*/
    }
    
                
        
    if ( _doTriplets){
        LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   Evaluating triplets" << flush; 
        // get triplet BSE Hamiltonian from _orbitalsAB
        ub::matrix<real_gwbse> _Hamiltonian_AB = _eh_d;
        
     
        
        
        //ub::matrix<real_gwbse>& _bseA= _orbitalsA->BSETripletCoefficients();
        
        
        const ub::matrix<real_gwbse>& _bseA = ub::project( _orbitalsA->BSETripletCoefficients(),
                ub::range (0, _orbitalsA->BSETripletCoefficients().size1() ), ub::range ( 0, _FeA )  );
        const ub::matrix<real_gwbse>& _bseB = ub::project( _orbitalsB->BSETripletCoefficients(),
                ub::range (0, _orbitalsB->BSETripletCoefficients().size1() ), ub::range ( 0, _FeB )  );
        //ub::matrix<real_gwbse>& _bseB = _orbitalsB->BSETripletCoefficients();
        
  
        
        
        bool _triplets = ProjectExcitons( _kap, _kbp,ctAB,ctBA, _bseA, _bseB, _Hamiltonian_AB, JAB_triplet);
        if ( _triplets ) {
            LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()   << "   calculated triplet couplings " << flush;
        }
    }
    
    
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Done with exciton couplings" << flush;
    return true;   
};


bool BSECoupling::ProjectExcitons(const ub::matrix<real_gwbse>& _kap,const ub::matrix<real_gwbse>& _kbp,
                                  const ub::matrix<real_gwbse>& ctAB,const ub::matrix<real_gwbse>& ctBA, 
                                  const ub::matrix<real_gwbse>& _bseA, const ub::matrix<real_gwbse>& _bseB, 
                                  const ub::matrix<real_gwbse>& _H, ub::matrix<real_gwbse>& _J){
    
    
    //cout << " Dimensions of _bseA " << _bseA.size1() << " : " << _bseA.size2() << endl;

     
     
     // get projection of monomer excitons on dimer product functions
     ub::matrix<real_gwbse> _proj_excA = ub::prod( ub::trans( _bseA ), _kap);
     ub::matrix<real_gwbse> _proj_excB = ub::prod( ub::trans( _bseB ), _kbp);
    
 
     //cout << "_proj_excA"<<_proj_excA.size1()<<"x"<<_proj_excA.size2()<<endl;
 
     //cout << "_proj_excB"<<_proj_excB.size1()<<"x"<<_proj_excB.size2()<<endl;
     
     
     //cout << "_ctAB"<<ctAB.size1()<<"x"<<ctAB.size2()<<endl;
     
     //cout << "_ctBA"<<ctBA.size1()<<"x"<<ctBA.size2()<<endl;
     
     unsigned _bseA_exc = _proj_excA.size1();
     unsigned _bseB_exc = _proj_excB.size1();
     unsigned _bse_exc=_bseA_exc+_bseB_exc;
     
     _J=ub::zero_matrix<real_gwbse>(_bse_exc,_bse_exc);
     unsigned _ctAB=ctAB.size1();
     
     unsigned _ctBA=ctBA.size1();
     unsigned _ct=_ctAB+_ctBA;
     unsigned nobasisfunc=_H.size1();
     
     
     
     ub::matrix<double> fe_states=ub::matrix<double>(_bse_exc,nobasisfunc);
     ub::project(fe_states, ub::range ( 0, _bseA_exc ),ub::range (0, nobasisfunc )  )=_proj_excA;
     ub::project(fe_states, ub::range ( _bseA_exc, _bse_exc ),ub::range (0, nobasisfunc )  )=_proj_excB;
      
     ub::matrix<double> ct_states=ub::matrix<double>(_ct,nobasisfunc);
     
     //cout<< _ct<< "ct states"<<endl;
    if(_ct>0){ 
     //orthogonalize ct-states with respect to the FE states. 
       LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << " Orthogonalizing CT-states with respect to FE-states" << flush;
   
     if(_ctAB>0){
     ub::project(ct_states, ub::range ( 0 , _ctAB ) ,ub::range (0,nobasisfunc ) )=ctAB;
    }
    if(_ctBA>0){
     ub::project(ct_states, ub::range ( _ctAB, _ct ),ub::range (0, nobasisfunc )  )=ctBA;
     }
       
       
        //orthogonalize ct-states with respect to FE states
     ub::matrix<double> overlaps=ub::prod(fe_states,ub::trans(ct_states));
     
     //cout << "overlap"<< overlaps.size1()<< "x"<<overlaps.size2()<<endl;
     //cout << overlaps<<endl;
     ub::matrix<double> correction=ub::prod(ub::trans(overlaps),fe_states);
     //cout << "correction"<< correction.size1()<< "x"<<correction.size2()<<endl;
    // cout << "ct_states"<< ct_states.size1()<< "x"<<ct_states.size2()<<endl;
  
   

         ct_states=ct_states-correction;    

    
    
     overlaps.resize(0,0);
     correction.resize(0,0);
     //normalize
    
     for (unsigned i=0;i<_ct;i++){
         double norm=0.0;
         for (unsigned j=0;j<nobasisfunc;j++){
         norm+=ct_states(i,j)*ct_states(i,j);    
         }
         //cout << "norm ["<<i<<"]:" <<norm<<endl;
         norm=1/std::sqrt(norm);
         if(norm<0.95){
            LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << " WARNING: CT-state "<< i<< " norm is only"<< norm << flush; 
         }
         for (unsigned j=0;j<nobasisfunc;j++){
            ct_states(i,j)=norm*ct_states(i,j);    
         }
         
     }
    //cout <<ub::prod(fe_states,ub::trans(ct_states))<<endl; 
      
    } 

     
     
     //cout << _ctAB <<_ctBA<<endl;
     //ub::matrix<double> _J_dimer = ub::zero_matrix<double>( _bse_exc +_ct, _bse_exc+_ct );
     //ub::matrix<double> _S_dimer = ub::zero_matrix<double>( _bse_exc +_ct, _bse_exc +_ct);
     
     //cout << _J_dimer.size1()<< " : "<<_J_dimer.size2()<<endl;
     //cout << _S_dimer.size1()<< " : "<<_S_dimer.size2()<<endl;
     
     #if (GWBSE_DOUBLE)
      const ub::matrix<double>& Htemp=_H;
#else
     ub::matrix<double> Htemp=_H;
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << " casting Hamiltonian to double  " << flush;

#endif
    
     ub::matrix<double> projection =ub::zero_matrix<double>(_bse_exc+_ct,nobasisfunc);
     
     
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << " merging projections into one vector  " << flush;
    
  ub::project(projection, ub::range (0 , _bse_exc) ,ub::range (0,nobasisfunc ) )=fe_states;
   
     if(_ct>0){
    ub::project(projection, ub::range ( _bse_exc , _bse_exc+_ct ) ,ub::range (0,nobasisfunc ) )=ct_states;
     }
      LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   Setting up coupling matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;
     // matrix _J
     
    //  E_A         J_AB        J_A_ABCT        J_A_BACT
    //  J_BA        E_B         J_B_ABCT        J_B_BACT
    //  J_ABCT_A    J_ABCT_B    E_ABCT          J_ABCT_BACT
    //  J_BACT_A   J_BACT_B    J_BACT_ABCT     E_BACT
     
     // I think this only works for hermitian/symmetric H so only in TDA
     // setup J
     
    
     
     ub::matrix<double> _temp=ub::prod(Htemp,ub::trans(projection));
     ub::matrix<double> _J_dimer=ub::prod(projection,_temp);
       #if (GWBSE_DOUBLE)
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << " no Casting  " << flush;
    #else
     Htemp.resize(0,0);
     #endif
     _temp.resize(0,0);
     

    
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   Setting up overlap matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;
     // setup S
    
    ub::matrix<double> _S_dimer=ub::prod(projection,ub::trans(projection));
    
    projection.resize(0,0);
    if(votca::tools::globals::verbose &&  _bse_exc+_ct<100){
         LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << "_J_dimer[Ryd]"<<flush;
     
     LOG(CTP::logDEBUG, *_pLog) << _J_dimer<<flush;
     LOG(CTP::logDEBUG, *_pLog) << "_S_dimer"<<flush;
     
     LOG(CTP::logDEBUG, *_pLog) << _S_dimer<<flush;
      LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    }
    //LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  Diagonlaize overlap "<< flush;
    //ub::vector<double> _S_eigenvalues;
     //cout << "_J_ortho"<<endl;
     
     //cout << _J_dimer<<endl;
     //linalg_eigenvalues(_S_eigenvalues,_S_dimer);
     //LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "  done "<< flush;
    double small=linalg_loewdin(_J_dimer,_S_dimer);
    
    if(votca::tools::globals::verbose && _bse_exc+_ct<100){
         LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    LOG(CTP::logDEBUG, *_pLog) << "_J_ortho[Ryd]"<<flush;
    LOG(CTP::logDEBUG, *_pLog) << _J_dimer<<flush;
    LOG(CTP::logDEBUG, *_pLog) << "_S-1/2"<<flush;
    LOG(CTP::logDEBUG, *_pLog) << _S_dimer<<flush;
     LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    }
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   Smallest value of dimer overlapmatrix is "<< small<< flush;
     
     //Diagonalize ct states
     
     if(_do_perturbation){
         bool _diag_ct=true;
     if (_ct > 0 && _diag_ct) {
        
         
                ub::matrix<double> transformation = ub::identity_matrix<double>(_bse_exc + _ct, _bse_exc + _ct);
                ub::vector<double> eigenvalues_ct;



                ub::matrix<double> Ct = ub::project(_J_dimer, ub::range(_bse_exc, _bse_exc + _ct), ub::range(_bse_exc, _bse_exc + _ct));
                linalg_eigenvalues(eigenvalues_ct, Ct);
                ub::project(transformation, ub::range(_bse_exc, _bse_exc + _ct), ub::range(_bse_exc, _bse_exc + _ct)) = Ct;
                
                Ct.resize(0, 0);

                if (votca::tools::globals::verbose) {

                    LOG(CTP::logDEBUG, *_pLog) << "FE state hamiltonian" << flush;
                    LOG(CTP::logDEBUG, *_pLog) << ub::project(_J_dimer, ub::range(0, _bse_exc), ub::range(0, _bse_exc)) << flush;
                    if (_ct > 0 ) {
                        LOG(CTP::logDEBUG, *_pLog) << "eigenvalues of CT states" << flush;
                        LOG(CTP::logDEBUG, *_pLog) << eigenvalues_ct << flush;
                    }

                }

                _temp = ub::prod(_J_dimer, transformation);
                _J_dimer = ub::prod(ub::trans(transformation), _temp);
                if (votca::tools::globals::verbose && _bse_exc + _ct < 100) {
                    LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------" << flush;
                    LOG(CTP::logDEBUG, *_pLog) << "_J_ortho[Ryd] CT-state diag" << flush;
                    LOG(CTP::logDEBUG, *_pLog) << _J_dimer << flush;
                    LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------" << flush;
                    }
     }
                    for (int stateA = 0; stateA < _levA; stateA++) {
                        double Ea = _J_dimer(stateA, stateA);
                        for (int stateB = 0; stateB < _levB; stateB++) {
                            int stateBd = stateB + _bseA_exc;
                            LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp() << "   Calculating coupling between exciton A" << stateA + 1 << " and exciton B" << stateB + 1 << flush;
                            double J = _J_dimer(stateA, stateBd);
                            
                            double Eb = _J_dimer(stateBd, stateBd);
                            for (unsigned k = _bse_exc; k < (_bse_exc + _ct); k++) {
                                double Eab = _J_dimer(k, k);
                                if (std::abs(Eab-Ea)<0.0001){
                                    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp() << "Energydifference between state A"<< stateA+1<< "and CT state "<< k+1<< " is "<< Eab-Ea<< "[Ryd]"<< flush;
                                }
                                if (std::abs(Eab-Eb)<0.0001){
                                    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp() << "Energydifference between state B"<< stateB+1<< "and CT state "<< k+1<< " is "<< Eab-Ea<< "[Ryd]"<< flush;

                                }
                                J += 0.5*_J_dimer(k, stateA) * _J_dimer(k, stateBd)*(1 / (Ea - Eab) + 1 / (Eb - Eab));// Have no clue why 0.5
                            }
                            _J(stateA, stateBd) = J;
                            _J(stateBd, stateA) = J;


                        }
                    }
     }
    
     
     
     else if(_do_full_diag){
     

     ub::vector<double> _J_eigenvalues;
   
  
     linalg_eigenvalues(_J_eigenvalues,_J_dimer);
     if(votca::tools::globals::verbose && _bse_exc+_ct<100){
          LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << "Eigenvectors of J"<<flush;
     
     LOG(CTP::logDEBUG, *_pLog) << _J_dimer<<flush;
     LOG(CTP::logDEBUG, *_pLog) << "J_eigenvalues[Ryd]"<< flush;
     LOG(CTP::logDEBUG, *_pLog) << _J_eigenvalues<< flush;
      LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     }
     //Calculate projection on subspace for every pair of excitons separately
     for (int stateA=0;stateA<_levA; stateA++){
          for (int stateB=0;stateB<_levB; stateB++){
              int stateBd=stateB+_bseA_exc;
              LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   Calculating coupling between exciton A"<< stateA+1<<" and exciton B"<<stateB+1 << flush;
              std::vector<unsigned> index;
              std::vector<int> signvec;
              for (unsigned i = 0; i < _bse_exc+_ct; i++) {
                        if (i == unsigned(stateA) || i == unsigned(stateBd)) {

                            double close = 0.0;
                            unsigned ind = 0;
                            int sign=0;
                            //row
                            for (unsigned j = 0; j < _bse_exc+_ct; j++) {
                                bool check = true;
                                // if index i is already in index
                                // should not happen but if one vector was similar to two others.
                                for (unsigned l = 0; l < index.size(); l++) {
                                    if (j == index[l]) {
                                        check = false;
                                        break;
                                    }
                                }

                                if (check && std::abs(_J_dimer(i, j)) > close) {
                                    ind = j;
                                    close = std::abs(_J_dimer(i, j));
                                    if (_J_dimer(i, j)>=0){
                                        sign=1;
                                    }
                                    else{
                                        sign=-1;
                                    }
                                }
                            }
                            index.push_back(ind);
                            signvec.push_back(sign);
                        }
                    }
              
              
    LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   Order is: [Initial state n->nth eigenvalue]"<<flush;
    LOG(CTP::logDEBUG, *_pLog) << "    A" << stateA+1<<":"<<stateA+1<<"->"<<index[0]+1<<" " ;   
    LOG(CTP::logDEBUG, *_pLog) << "    B" << stateB+1<<":"<<stateBd+1<<"->"<<index[1]+1<<" "<<flush ;  
       
   
     
     //setting up transformation matrix _T and diagonal matrix _E for the eigenvalues;
     
     ub::matrix<double> _E=ub::zero_matrix<double>(2,2);
     ub::matrix<double> _T=ub::zero_matrix<double>(2,2);
     //find the eigenvectors which are most similar to the initial states
     
    //row 
    for (unsigned i = 0; i < 2; i++) {
    unsigned k=index[i];
    double sign=signvec[i];       
    double normr=1/std::sqrt(_J_dimer(stateA,k)*_J_dimer(stateA,k)+_J_dimer(stateBd,k)*_J_dimer(stateBd,k));
           _T(0,i ) = sign*_J_dimer(stateA,k)* normr;
           _T(1,i ) = sign*_J_dimer(stateBd,k)*normr;
           _E(i,i)=_J_eigenvalues(k);
       }
     
   
     if ((_T(1,1)*_T(0,0)-_T(1,0)*_T(0,1))<0){
         LOG(CTP::logDEBUG, *_pLog) << " Reduced state matrix is not in a right handed basis, multiplying second eigenvector by -1 "<<flush ;
         _T(0,1)=-_T(0,1);
         _T(1,1)=-_T(1,1);
     }
     
      if(votca::tools::globals::verbose){
     LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << "_T"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << _T<<flush;
     
     }
     
     ub::matrix<double> S_small=ub::prod(_T,ub::trans(_T));
      if(votca::tools::globals::verbose){
     
     LOG(CTP::logDEBUG, *_pLog) << "S_small"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << S_small<<flush;
    
     }
     //orthogonalize that matrix
     small=linalg_loewdin(_E,S_small);
     LOG(CTP::logDEBUG, *_pLog) << CTP::TimeStamp()  << "   Smallest value of dimer overlapmatrix is "<< small<< flush;
       if(votca::tools::globals::verbose){
    
     LOG(CTP::logDEBUG, *_pLog) << "S-1/2"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << S_small<<flush;
      LOG(CTP::logDEBUG, *_pLog) << "E_ortho"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << _E<<flush;
     }
     _T=ub::prod(_T,S_small);
     //cout <<  ub::prod(_T,ub::trans(_T))<<endl;
   if(votca::tools::globals::verbose){
    
     LOG(CTP::logDEBUG, *_pLog) << "T_ortho"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << _T<<flush;
    LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     }
     
     ub::matrix<double> temp=ub::prod(_T,_E);
     
      
    
    
     ub::matrix<double> _J_small=ub::prod(temp,ub::trans(_T));
     if(votca::tools::globals::verbose){
    
     LOG(CTP::logDEBUG, *_pLog) << "T_ortho*E_ortho"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << temp<<flush;
      LOG(CTP::logDEBUG, *_pLog) << "T_ortho*E_ortho*T_ortho^T"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << _J_small<<flush;
     }
     
          
     _J(stateA,stateBd)=_J_small(0,1);
     //_J(stateA,stateBd)=_J_small(0,0);
    // _J(stateBd,stateBd)=_J_small(1,1);
     _J(stateBd,stateA)=_J_small(1,0);
     //cout << _temp2<<endl;
   
     
    
         }
     }
  
     }
       if(votca::tools::globals::verbose){
     LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << "Jeff[Ryd]"<<flush;
     LOG(CTP::logDEBUG, *_pLog) << _J<<flush;
     LOG(CTP::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     }
      
     return true;
    
    
    
    
}
    
}}
