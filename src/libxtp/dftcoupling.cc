/*
 *            Copyright 2009-2018 The VOTCA Development Team
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


#include <votca/xtp/aomatrix.h>
#include <votca/xtp/dftcoupling.h>

#include <votca/tools/constants.h>
#include <boost/format.hpp>
#include <boost/progress.hpp>


namespace votca { namespace xtp {


  using std::flush;
using boost::format;


double DFTcoupling::getCouplingElement( int levelA, int levelB,  Orbitals* _orbitalsA,
    Orbitals* _orbitalsB, Eigen::MatrixXd* _JAB, double  _energy_difference ) {

    
    int _levelsA = _orbitalsA->getNumberOfLevels();
    
    if ( _energy_difference != 0 ) {
        std::vector<int> list_levelsA = *_orbitalsA->getDegeneracy( levelA, _energy_difference );
        std::vector<int> list_levelsB = *_orbitalsA->getDegeneracy( levelB, _energy_difference );
        
        double _JAB_sq = 0; double _JAB_one_level;
        
        for (std::vector<int>::iterator iA = list_levelsA.begin()++; iA != list_levelsA.end(); iA++) {
                for (std::vector<int>::iterator iB = list_levelsB.begin()++; iB != list_levelsB.end(); iB++) { 
                    _JAB_one_level = (*_JAB)( *iA - 1  , *iB -1 + _levelsA );
                    _JAB_sq +=  _JAB_one_level*_JAB_one_level ;
                }
        }
        
        return sqrt(_JAB_sq / ( list_levelsA.size() * list_levelsB.size() ) ) * tools::conv::hrt2ev ;
        
    } else {
        
        return (*_JAB)( levelA - 1  , levelB -1 + _levelsA ) * tools::conv::hrt2ev;
        
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
    Orbitals* _orbitalsAB, Eigen::MatrixXd* _JAB) {

    CTP_LOG(ctp::logDEBUG,*_pLog) << "Calculating electronic couplings" << flush;
    
    const std::vector<QMAtom*> atomsA=_orbitalsA->QMAtoms();
    const std::vector<QMAtom*> atomsB=_orbitalsB->QMAtoms();
     const std::vector<QMAtom*> atomsAll = _orbitalsAB->QMAtoms();
     
  for (unsigned i = 0; i < atomsAll.size(); i++) {
      QMAtom* dimer = atomsAll[i];
      QMAtom* monomer = NULL;

      if (i < atomsA.size()) {
          monomer = atomsA[i];
      } else if (i < atomsB.size() + atomsA.size()) {
          monomer = atomsB[i - atomsA.size()]; 
      } else {
          // Linker
          CTP_LOG(ctp::logERROR, *_pLog) << (format("Neither Monomer A nor Monomer B contains atom %s on line %u. Hence, this atom is part of a linker.") %dimer->getType() %(i+1) ).str()<<flush;
          continue;
      }
      
      if(!monomer->getPos().isClose(dimer->getPos(), 0.001)){
              CTP_LOG(ctp::logINFO, *_pLog) << "======WARNING=======\n Coordinates of monomer and dimer atoms do not agree, do you know what you are doing?" << flush;
          }

      if (monomer->getType() != dimer->getType()) {
          throw std::runtime_error("\nERROR: Atom types do not agree in dimer and monomers\n");
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
    //       | 0          Orbitals_B |.T  X   | Overlap_B |  X  ( Orbitals_AB )
 
    
    Eigen::MatrixXd _psi_AxB=Eigen::MatrixXd::Zero( _basisA + _basisB, _levelsA + _levelsB  );
    
      CTP_LOG(ctp::logDEBUG,*_pLog) << "Constructing direct product AxB [" 
            << _psi_AxB.rows() << "x" 
            << _psi_AxB.cols() << "]" << flush;    
    
    // constructing merged orbitals
    _psi_AxB.block(0,0, _basisA,_levelsA) = _orbitalsA->MOCoefficients();
    _psi_AxB.block(_basisA,_levelsA, _basisB,_levelsB) =_orbitalsB->MOCoefficients(); 
   
    Eigen::MatrixXd overlap;
    if ( _orbitalsAB->hasAOOverlap() ) {
            CTP_LOG(ctp::logDEBUG,*_pLog) << "Reading overlap matrix from orbitals" << flush; 
           overlap= _orbitalsAB->AOOverlap();
    }else{
        CTP_LOG(ctp::logDEBUG,*_pLog) << "Calculating overlap matrix for basisset: "<< _orbitalsAB->getDFTbasis()<< flush; 
        BasisSet _dftbasisset;
        AOBasis _dftbasis;
        _dftbasisset.LoadBasisSet(_orbitalsAB->getDFTbasis());

        _dftbasis.AOBasisFill(_dftbasisset, _orbitalsAB->QMAtoms());
        AOOverlap _dftAOoverlap;
        _dftAOoverlap.Fill(_dftbasis);
        overlap=_dftAOoverlap.Matrix();
    }
     CTP_LOG(ctp::logDEBUG,*_pLog) << "Projecting dimer onto monomer orbitals" << flush; 
    Eigen::MatrixXd _psi_AxB_dimer_basis =_psi_AxB.transpose()*overlap*_orbitalsAB->MOCoefficients(); 

    unsigned int LevelsA = _levelsA;
    for (unsigned i=0;i<_psi_AxB_dimer_basis.rows();i++){
        double mag=_psi_AxB_dimer_basis.row(i).squaredNorm();
        if (mag<0.95){
            int monomer = 0;
            int level = 0;
            if ( i < LevelsA ) {
                monomer = 1;
                level   = i;
            } else {
                monomer = 2;
                level   = i -_levelsA;
                
            }
            CTP_LOG(ctp::logERROR,*_pLog) << "\nERROR: " << i << " Projection of orbital " << level << " of monomer " << monomer << " on dimer is insufficient,mag="<<mag<<" maybe the orbital order is screwed up, otherwise increase dimer basis.\n"<<flush;
        }
    }
    // J = psi_AxB_dimer_basis * FAB * psi_AxB_dimer_basis^T
    CTP_LOG(ctp::logDEBUG,*_pLog) << "Projecting the Fock matrix onto the dimer basis" << flush;   
      
    Eigen::MatrixXd JAB_dimer = _psi_AxB_dimer_basis*_orbitalsAB->MOEnergies().asDiagonal()*_psi_AxB_dimer_basis.transpose();  
    // S = psi_AxB_dimer_basis * psi_AxB_dimer_basis^T
    CTP_LOG(ctp::logDEBUG,*_pLog) << "Constructing Overlap matrix" << flush;    
    Eigen::MatrixXd _S_AxB = _psi_AxB_dimer_basis*_psi_AxB_dimer_basis.transpose();  
     
   CTP_LOG(ctp::logDEBUG,*_pLog) << "Calculating the effective overlap JAB [" 
              << JAB_dimer.rows() << "x" 
              << JAB_dimer.cols() << "]" << flush;  
   
   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_S_AxB);
   Eigen::MatrixXd Sm1=es.operatorInverseSqrt();
   CTP_LOG(ctp::logDEBUG,*_pLog) << "Smallest eigenvalue of overlap matrix is "<<es.eigenvalues()(0)<< flush;    
   (*_JAB) = Sm1*JAB_dimer*Sm1;
    
    CTP_LOG(ctp::logDEBUG,*_pLog) << "Done with electronic couplings" << flush;
    
    return true;   

}


    
}}
