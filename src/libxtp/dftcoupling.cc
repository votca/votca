/*
 *            Copyright 2009-2017 The VOTCA Development Team
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
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/dftcoupling.h>

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/constants.h>
#include <boost/format.hpp>
#include <boost/progress.hpp>


namespace votca { namespace xtp {

namespace ub = boost::numeric::ublas;


using boost::format;


double DFTcoupling::getCouplingElement( int levelA, int levelB,  Orbitals* _orbitalsA,
    Orbitals* _orbitalsB, ub::matrix<double>* _JAB, double  _energy_difference ) {

    
    int _levelsA = _orbitalsA->getNumberOfLevels();
    
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
        
        return sqrt(_JAB_sq / ( list_levelsA.size() * list_levelsB.size() ) ) * tools::conv::hrt2ev ;
        
    } else {
        
        return _JAB->at_element( levelA - 1  , levelB -1 + _levelsA ) * tools::conv::hrt2ev;
        
    }
    // the  matrix should be symmetric, could also return this element
    // _JAB.at_element( _levelsA + levelB - 1  , levelA - 1 );
}

/**
 * \brief evaluates electronic couplings  
 * 
 * This is a fast version with a rather large block matrix AxB
 * 
 * @param _orbitalsA molecular orbitals of molecule A
 * @param _orbitalsB molecular orbitals of molecule B
 * @param _orbitalsAB molecular orbitals of the dimer AB
 * @param _JAB matrix with electronic couplings
 * @return false if failed
 */
bool DFTcoupling::CalculateIntegrals(Orbitals* _orbitalsA, Orbitals* _orbitalsB, 
    Orbitals* _orbitalsAB, ub::matrix<double>* _JAB) {

    CTP_LOG(ctp::logDEBUG,*_pLog) << "Calculating electronic couplings" << flush;
    
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
            throw runtime_error((format("Number of Atoms in dimer %3i and the two monomers A:%3i B:%3i does not agree") %atomsAB.size() %atomsA.size() %atomsB.size()).str());
        }
        
      if(monomer->getType() != dimer->getType()){
            throw runtime_error("\nERROR: Atom types do not agree in dimer and monomers\n");
        }
        if(tools::abs(monomer->getPos()-dimer->getPos())>0.001){
            CTP_LOG(ctp::logERROR,*_pLog) << "======WARNING=======\n Coordinates of monomers and dimer atoms do not agree, do you know what you are doing?\n " << flush;
            break;
        }
        
    }
    
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        CTP_LOG(ctp::logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
        return false;
    }
        
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    //boost::timer t; // start timing
    //double _st = t.elapsed();
    
    CTP_LOG(ctp::logDEBUG,*_pLog) << "Levels:Basis A[" << _levelsA << ":" << _basisA << "]"
                                     << " B[" << _levelsB << ":" << _basisB << "]" << flush;
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        CTP_LOG(ctp::logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
        return false;
    } 
    
    //       | Orbitals_A          0 |      | Overlap_A |     
    //       | 0          Orbitals_B |  X   | Overlap_B |  X  Transpose( Orbitals_AB )
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
    ub::matrix<double> _psi_AxB ( _levelsA + _levelsB, _basisA + _basisB  );
    

    CTP_LOG(ctp::logDEBUG,*_pLog) << "Constructing direct product AxB [" 
            << _psi_AxB.size1() << "x" 
            << _psi_AxB.size2() << "]" << flush;    
    
    
    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = _orbitalsA->MOCoefficients();
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = _orbitalsB->MOCoefficients(); 

    // psi_AxB * S_AB * psi_AB
    
    
    CTP_LOG(ctp::logDEBUG,*_pLog) << "Projecting dimer onto monomer orbitals" << flush; 
    ub::matrix<double> overlap;
    if ( _orbitalsAB->hasAOOverlap() ) {
            CTP_LOG(ctp::logDEBUG,*_pLog) << "Reading overlap matrix from orbitals" << flush; 
           overlap= _orbitalsAB->AOOverlap();
    }else{
        CTP_LOG(ctp::logDEBUG,*_pLog) << "Calculating overlap matrix for basisset: "<< _orbitalsAB->getDFTbasis()<< flush; 
        BasisSet _dftbasisset;
        AOBasis _dftbasis;
        _dftbasisset.LoadBasisSet(_orbitalsAB->getDFTbasis());

        _dftbasis.AOBasisFill(&_dftbasisset, _orbitalsAB->QMAtoms());
        AOOverlap _dftAOoverlap;
        _dftAOoverlap.Fill(_dftbasis);
        overlap=_dftAOoverlap.Matrix();
    }
    
    ub::matrix<double> _psi_AB = ub::prod(overlap, ub::trans( _orbitalsAB->MOCoefficients()) ); 
    ub::matrix<double> _psi_AxB_dimer_basis = ub::prod( _psi_AxB, _psi_AB );  
    _psi_AB.clear();
    
    //check to see if projection quality is sufficient
    for (unsigned i=0;i<_psi_AxB_dimer_basis.size1();i++){
        double mag=0.0;
        for (unsigned j=0;j<_psi_AxB_dimer_basis.size2();j++){
            mag+=_psi_AxB_dimer_basis(i,j)*_psi_AxB_dimer_basis(i,j);
            
    }
        if (mag<0.95){
            throw runtime_error("\nERROR: Projection of monomer orbitals on dimer is insufficient, maybe the orbital order is screwed up, otherwise increase dimer basis.\n");
        }
    }
 
    // J = psi_AxB_dimer_basis * FAB * psi_AxB_dimer_basis^T
    CTP_LOG(ctp::logDEBUG,*_pLog) << "Projecting the Fock matrix onto the dimer basis" << flush;   
    ub::diagonal_matrix<double> _fock_AB( _orbitalsAB->getNumberOfLevels(), _orbitalsAB->MOEnergies().data() ); 
    ub::matrix<double> _temp = ub::prod( _fock_AB, ub::trans( _psi_AxB_dimer_basis ) ) ; 
    ub::matrix<double> JAB_dimer = ub::prod( _psi_AxB_dimer_basis, _temp);  
 
    // S = psi_AxB_dimer_basis * psi_AxB_dimer_basis^T
    CTP_LOG(ctp::logDEBUG,*_pLog) << "Constructing Overlap matrix" << flush;    
    ub::matrix<double> _S_AxB = ub::prod( _psi_AxB_dimer_basis, ub::trans( _psi_AxB_dimer_basis ));  
     
   CTP_LOG(ctp::logDEBUG,*_pLog) << "Calculating the effective overlap JAB [" 
              << JAB_dimer.size1() << "x" 
              << JAB_dimer.size2() << "]" << flush;  
   double smalleig=linalg_loewdin(JAB_dimer, _S_AxB);
    CTP_LOG(ctp::logDEBUG,*_pLog) << "Smallest eigenvalue of overlap matrix is "<<smalleig<< flush;    
    (*_JAB) = JAB_dimer;    
    
    
    CTP_LOG(ctp::logDEBUG,*_pLog) << "Done with electronic couplings" << flush;
    
    return true;   

}


    
}}
