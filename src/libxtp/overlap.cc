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

#include <votca/xtp/overlap.h>
#include <votca/tools/linalg.h>

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

double inv_sqrt(double x) { return 1./sqrt(x); }
using boost::format;

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
        
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    //boost::timer t; // start timing
    //double _st = t.elapsed();
    
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
            << _psi_AxB.size2() << "]" << flush;    
    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = *_orbitalsA->getOrbitals();
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = *_orbitalsB->getOrbitals(); 

    // psi_AxB * S_AB * psi_AB
    LOG(logDEBUG,*_pLog) << "Projecting dimer onto monomer orbitals" << flush; 
    ub::matrix<double> _orbitalsAB_Transposed = ub::trans( *_orbitalsAB->getOrbitals() );  
    if ( (*_orbitalsAB->getOverlap()).size1() == 0 ) {
            LOG(logERROR,*_pLog) << "Overlap matrix is not stored"; 
            return false;
    }
     
    ub::matrix<double> _psi_AB = ub::prod( *_orbitalsAB->getOverlap(), _orbitalsAB_Transposed );  
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
    LOG(logDEBUG,*_pLog) << "Projecting the Fock matrix onto the dimer basis" << flush;   
    ub::diagonal_matrix<double> _fock_AB( _orbitalsAB->getNumberOfLevels(), (*_orbitalsAB->getEnergies()).data() ); 
    ub::matrix<double> _temp = ub::prod( _fock_AB, ub::trans( _psi_AxB_dimer_basis ) ) ; 
    ub::matrix<double> JAB_dimer = ub::prod( _psi_AxB_dimer_basis, _temp);  
 
    // S = psi_AxB_dimer_basis * psi_AxB_dimer_basis^T
    LOG(logDEBUG,*_pLog) << "Constructing Overlap matrix" << flush;    
    ub::symmetric_matrix<double> _S_AxB = ub::prod( _psi_AxB_dimer_basis, ub::trans( _psi_AxB_dimer_basis ));  
    ub::matrix<double> _S_AxB_2(_S_AxB.size1(), _S_AxB.size1() );
    ub::trans( _S_AxB );
  
    // Square root of the overlap matrix
    LOG(logDEBUG,*_pLog) << "Calculating square root of the overlap matrix" << flush;    
    SQRTOverlap( _S_AxB , _S_AxB_2 );
 
     
   LOG(logDEBUG,*_pLog) << "Calculating the effective overlap JAB [" 
              << JAB_dimer.size1() << "x" 
              << JAB_dimer.size2() << "]" << flush;  
       
    ub::matrix<double> JAB_temp( _levelsA + _levelsB, _levelsA + _levelsB ); 
    ub::noalias(JAB_temp) = ub::prod( JAB_dimer, _S_AxB_2 );  
    (*_JAB) = ub::prod( _S_AxB_2, JAB_temp );    
    
    
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
