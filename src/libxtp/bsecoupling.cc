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
#include <votca/xtp/bsecoupling.h>
#include <votca/tools/constants.h>
#include <boost/format.hpp>



namespace votca { namespace xtp {

using boost::format;
using namespace tools;

void BSECoupling::Initialize(Property* options){
    
    #if (GWBSE_DOUBLE)
        XTP_LOG(xtp::logDEBUG, *_pLog) <<  " Compiled with full double support" << flush;   
    #else
        XTP_LOG(xtp::logDEBUG, *_pLog) <<  " Compiled with float/double mixture (standard)" << flush;   
    #endif
    
    std::string key = Identify(); 
    _doSinglets=false;
    _doTriplets=false;
   _output_perturbation=false;
    
    
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
     
     
   if ( options->exists( key + ".algorithm") ) {
                string algorithm = options->get(key + ".algorithm").as<string> ();
                 if(algorithm=="perturbation"){
                    _output_perturbation=true;
                 }
   }
    
   
    
        _levA  = options->get(key + ".moleculeA.states").as<int> ();
        _levB  = options->get(key + ".moleculeB.states").as<int> ();
        _occA  = options->get(key + ".moleculeA.occLevels").as<int> ();
        _occB  = options->get(key + ".moleculeB.occLevels").as<int> ();
        _unoccA  = options->get(key + ".moleculeA.unoccLevels").as<int> ();
        _unoccB  = options->get(key + ".moleculeB.unoccLevels").as<int> ();
        
   
        
        
}

void BSECoupling::addoutput(Property *_type_summary,Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB){
   
    string algorithm="full_diag";
    int methodindex=1;
    if (_output_perturbation){
        algorithm="perturbation";
        methodindex=0;
    }
    _type_summary->setAttribute("algorithm",algorithm);
    if (_doSinglets){
        Property *_singlet_summary = &_type_summary->add("singlets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB <_levB ; ++stateB ) {
               double JAB = getSingletCouplingElement( stateA , stateB, methodindex);
              
               Property *_coupling_summary = &_singlet_summary->add("coupling", (format("%1$1.6e") % JAB).str()); 
               double energyA = _orbitalsA->BSESingletEnergies()(stateA)*conv::hrt2ev;
               double energyB = _orbitalsB->BSESingletEnergies()(stateB)*conv::hrt2ev;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", (format("%1$1.6e") % energyA).str());
               _coupling_summary->setAttribute("energyB", (format("%1$1.6e") % energyB).str());
               _coupling_summary->setAttribute("pert", (format("%1$1.6e") % getSingletCouplingElement( stateA , stateB, 0)).str());
               _coupling_summary->setAttribute("diag", (format("%1$1.6e") % getSingletCouplingElement( stateA , stateB, 1)).str());
               
           } 
        }
    }
    
    
    if ( _doTriplets){
        Property *_triplet_summary = &_type_summary->add("triplets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB < _levA ; ++stateB ) {
               double JAB = getTripletCouplingElement( stateA , stateB,methodindex );
               //real_gwbse energyAD = getTripletDimerEnergy( stateA  );
               //real_gwbse energyBD = getTripletDimerEnergy( stateB  );
               Property *_coupling_summary = &_triplet_summary->add("coupling", (format("%1$1.6e") % JAB).str()); 
               double energyA = _orbitalsA->BSETripletEnergies()(stateA)*conv::hrt2ev;
               double energyB = _orbitalsB->BSETripletEnergies()(stateB)*conv::hrt2ev;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", (format("%1$1.6e") % energyA).str());
               _coupling_summary->setAttribute("energyB", (format("%1$1.6e") % energyB).str());
               _coupling_summary->setAttribute("pert", (format("%1$1.6e") % getTripletCouplingElement( stateA , stateB, 0)).str());
               _coupling_summary->setAttribute("diag", (format("%1$1.6e") % getTripletCouplingElement( stateA , stateB, 1)).str());
              
           } 
        }
    }       
}


double BSECoupling::getSingletCouplingElement( int levelA, int levelB, int methodindex) {
    return JAB_singlet[methodindex]( levelA  , levelB +  _levA ) * votca::tools::conv::hrt2ev;
}



double BSECoupling::getTripletCouplingElement( int levelA, int levelB, int methodindex) {

    return JAB_triplet[methodindex]( levelA  , levelB + _levA ) * votca::tools::conv::hrt2ev;
}


/**
 * \brief evaluates electronic couplings  
 *   
 * @param _orbitalsA molecular orbitals of molecule A
 * @param _orbitalsB molecular orbitals of molecule B
 * @param _orbitalsAB molecular orbitals of the dimer AB
 * @param _JAB matrix with electronic couplings
 * @return false if failed
 */
bool BSECoupling::CalculateCouplings(Orbitals* _orbitalsA, Orbitals* _orbitalsB, Orbitals* _orbitalsAB) {
       XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Calculating exciton couplings" << flush;
     // set the parallelization 
    #ifdef _OPENMP
    
    if ( _openmp_threads > 0 ) omp_set_num_threads(_openmp_threads);      
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << " Using "<< omp_get_max_threads()<<" threads" << flush;
    #endif
    
    
    
    _orbitalsAB->setCoupledExcitonsA(_levA);
    _orbitalsAB->setCoupledExcitonsB(_levB);
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
        
        if(monomer->getType() != dimer->getType()){
            throw runtime_error("\nERROR: Atom types do not agree in dimer and monomers\n");
        }
        if(tools::abs(monomer->getPos()-dimer->getPos())>0.001){
            XTP_LOG(xtp::logERROR,*_pLog) << "======WARNING=======\n Coordinates of monomers and dimer atoms do not agree, do you know what you are doing?\n " << flush;
            break;
        }
        
    }
    
    
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        XTP_LOG(xtp::logERROR,*_pLog) << "Basis set size is not stored in monomers" << flush;
        return false;
    }

    // number of levels stored in monomers
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
        
    // get exciton information of molecule A
    int _bseA_cmax        = _orbitalsA->getBSEcmax();
    int _bseA_cmin        = _orbitalsA->getBSEcmin();
    int _bseA_vmax        = _orbitalsA->getBSEvmax();
    int _bseA_vmin        = _orbitalsA->getBSEvmin();
    int _bseA_vtotal      = _bseA_vmax - _bseA_vmin +1 ;
    int _bseA_ctotal      = _bseA_cmax - _bseA_cmin +1 ;
    int _bseA_size        = _bseA_vtotal * _bseA_ctotal;
    int _bseA_singlet_exc = _orbitalsA->BSESingletCoefficients().cols();
    int _bseA_triplet_exc = _orbitalsA->BSETripletCoefficients().cols();

    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()   << "   molecule A has " << _bseA_singlet_exc << " singlet excitons with dimension " << _bseA_size << flush;
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()   << "   molecule A has " << _bseA_triplet_exc << " triplet excitons with dimension " << _bseA_size << flush;
    
    // now, two storage assignment matrices for two-particle functions
    Eigen::MatrixXi _combA;
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
    int _bseB_singlet_exc = _orbitalsB->BSESingletCoefficients().cols();
    int _bseB_triplet_exc = _orbitalsB->BSETripletCoefficients().cols();
    

    
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()   << "   molecule B has " << _bseB_singlet_exc << " singlet excitons with dimension " << _bseB_size << flush;
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()   << "   molecule B has " << _bseB_triplet_exc << " triplet excitons with dimension " << _bseB_size << flush;
    
    // now, two storage assignment matrices for two-particle functions
    Eigen::MatrixXi _combB;
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
        XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Number of excitons you want is greater than stored for molecule A. Setting to max number available" << flush; 
        _levA=_bseA_singlet_exc;
    }
    if(_levB>_bseB_singlet_exc){
        XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Number of excitons you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levB=_bseB_singlet_exc;
    }
    
    
    if(_levA>_bseA_singlet_exc){
        XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Number of Frenkel states you want is greater than stored for molecule A. Setting to max number available" << flush; 
        _levA=_bseA_singlet_exc;
    }
    if(_levB>_bseB_singlet_exc){
        XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Number of Frenkel states you want is greater than stored for molecule B. Setting to max number available" << flush; 
        _levB=_bseB_singlet_exc;
    }
    
    if(_unoccA>_bseA_ctotal){
        XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Number of occupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccA=_bseA_ctotal;
    }
    else if (_unoccA<0){
        _unoccA=_bseA_ctotal;
    }
    if(_unoccB>_bseB_ctotal){
        XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Number of occupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _unoccB=_bseB_ctotal;
    }
    else if (_unoccB<0){
        _unoccB=_bseB_ctotal;
    }
    
    if(_occA>_bseA_vtotal){
        XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Number of unoccupied orbitals in molecule A for CT creation exceeds number of KS-orbitals in BSE" << flush; 
        _occA=_bseA_vtotal;
    }
    else if (_occA<0){
        _occA=_bseA_vtotal;
    }
    if(_occB>_bseB_vtotal){
        XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Number of unoccupied orbitals in molecule B for CT creation exceeds number of KS-orbitals in BSE" << flush; 
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
    if ( ! _orbitalsAB->hasEHinteraction_triplet() && _doTriplets){
        XTP_LOG(xtp::logERROR,*_pLog) << "BSE EH for triplets not stored " << flush;
        return false;
    }
    if ( ! _orbitalsAB->hasEHinteraction_singlet() && _doSinglets){
        XTP_LOG(xtp::logERROR,*_pLog) << "BSE EH for singlets not stored " << flush;
        return false;
    }
    const MatrixXfd&    _eh_t = _orbitalsAB->eh_t(); 
    const MatrixXfd&    _eh_s = _orbitalsAB->eh_s(); 
    if(_doTriplets){
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()   << "   dimer AB has BSE EH interaction triplet with dimension " << _eh_t.rows() << " x " <<  _eh_t.cols() << flush;
    }
    if(_doSinglets){
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()   << "   dimer AB has BSE EH interaction singlet with dimension " << _eh_s.rows() << " x " <<  _eh_s.cols() << flush;
    }
    // now, two storage assignment matrices for two-particle functions
    Eigen::MatrixXi _combAB;
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
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "   levels used in BSE of molA: " << _bseA_vmin << " to " << _bseA_cmax << " total: " << _bseA_vtotal + _bseA_ctotal <<  flush;
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "   levels used in BSE of molB: " << _bseB_vmin << " to " << _bseB_cmax << " total: " << _bseB_vtotal + _bseB_ctotal <<  flush;
    
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        XTP_LOG(xtp::logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
        return false;
    } 
    
    //       | Orbitals_A          0 |      | Overlap_A |     
    //       | 0          Orbitals_B |.T  X   | Overlap_B |  X  ( Orbitals_AB )
    
    Eigen::MatrixXd _psi_AxB =Eigen::MatrixXd::Zero( _levelsA + _levelsB, _basisA + _basisB  );
    
  
    
    // constructing merged orbitals
    _psi_AxB.block(0,0,_levelsA , _basisA) = _orbitalsA->MOCoefficients().block(_bseA_vmin,0, _bseA_cmax+1-_bseA_vmin, _basisA );
    _psi_AxB.block(_levelsA, _basisA,_levelsB,_basisB) =_orbitalsB->MOCoefficients().block(_bseB_vmin,0,_bseB_cmax+1-_bseA_vmin,_basisB); 
    
    // psi_AxB * S_AB * psi_AB
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "   projecting monomer onto dimer orbitals" << flush; 
    
     Eigen::MatrixXd _overlapAB;
    if ( _orbitalsAB->hasAOOverlap() ) {
            XTP_LOG(xtp::logDEBUG,*_pLog) << "Reading overlap matrix from orbitals" << flush; 
           _overlapAB= _orbitalsAB->AOOverlap();
    }else{
        XTP_LOG(xtp::logDEBUG,*_pLog) << "Calculating overlap matrix for basisset: "<< _orbitalsAB->getDFTbasis()<< flush; 
        BasisSet _dftbasisset;
        AOBasis _dftbasis;
        _dftbasisset.LoadBasisSet(_orbitalsAB->getDFTbasis());

        _dftbasis.AOBasisFill(&_dftbasisset, _orbitalsAB->QMAtoms());
        AOOverlap _dftAOoverlap;
        _dftAOoverlap.Fill(_dftbasis);
        _overlapAB=_dftAOoverlap.Matrix();
    }
    
  
    Eigen::MatrixXd _psi_AxB_dimer_basis = _psi_AxB.transpose()*_overlapAB*_orbitalsAB->MOCoefficients();  
    _overlapAB.resize(0,0);
    
    
    //cout<< "_psi_AxB_dimer"<<endl;
    unsigned int LevelsA = _levelsA;
    for (unsigned i=0;i<_psi_AxB_dimer_basis.rows();i++){
        double mag=_psi_AxB_dimer_basis.row(i).squaredNorm();
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
            XTP_LOG(xtp::logERROR,*_pLog) << "\nERROR: " << i << " Projection of orbital " << level << " of monomer " << monomer << " on dimer is insufficient,mag="<<mag<<" maybe the orbital order is screwed up, otherwise increase dimer basis.\n"<<flush;
        }
    }
   
    
    //notation AB is CT states with A+B-, BA is the counterpart
    //Setting up CT-states:
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "   Setting up CT-states" << flush; 
    //Number of A+B- states
    int noAB=_occA*_unoccB;
    //Number of A-B+ states
    int noBA=_unoccA*_occB;
    
    
    Eigen::MatrixXi comb_CTAB=Eigen::MatrixXi::Zero(noAB,2);
    
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
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  <<"  "<<noAB <<" CT states A+B- created" << flush;
 
    Eigen::MatrixXi  comb_CTBA=Eigen::MatrixXi::Zero(noBA,2);
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
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  <<"  "<<noBA <<" CT states B+A- created" << flush;
    
    
    
    // these 4 matrixes, matrix(i,j) contains the j-th dimer MO component of the i-th excitation
   
    ctAB.resize(noAB,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_CT = 0 ; _i_CT < noAB ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctAB(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTAB(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTAB(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
    
    ctBA.resize(noBA,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_CT = 0 ; _i_CT < noBA ; _i_CT++){
    for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
        ctBA(_i_CT,_i_bseAB)=_psi_AxB_dimer_basis( comb_CTBA(_i_CT,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( comb_CTBA(_i_CT,1), _combAB( _i_bseAB,1) );
        }
    }
    
      
    // some more convenient storage
    
    
    _kap.resize(_bseA_size,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) );
            
        }
    }

    
    

    _kbp.resize(_bseB_size,_bseAB_size);
    #pragma omp parallel for
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) );
        }
    }
    
    // Same routines but also take <v|c'> <c|v'> projections into account 
    /*
 
    _kap.resize(_bseA_size,_bseAB_size);
    for ( int _i_bseA = 0 ; _i_bseA < _bseA_size ; _i_bseA++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kap(_i_bseA,_i_bseAB) = _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,1) )+
              _psi_AxB_dimer_basis( _combA(_i_bseA,1), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combA(_i_bseA,0), _combAB( _i_bseAB,1) )     ;
            
        }
    }

    

    
   
    _kbp.resize(_bseB_size,_bseAB_size);
    for ( int _i_bseB = 0 ; _i_bseB < _bseB_size ; _i_bseB++){
        for ( int _i_bseAB = 0 ; _i_bseAB < _bseAB_size ; _i_bseAB++){
            _kbp(_i_bseB,_i_bseAB) = _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,1) )+
                    _psi_AxB_dimer_basis( _combB(_i_bseB,1), _combAB( _i_bseAB,0) ) * _psi_AxB_dimer_basis( _combB(_i_bseB,0), _combAB( _i_bseAB,1) );
        }
    }
    */ 
  

    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()   << "   construct projection of product functions " << flush; 

   
    
        
    //     cout << "Size of _kap " << _kap.size1() << " : " <<  _kap.size2() << "\n" << flush; 
    //     cout << "Size of _kbp " << _kbp.size1() << " : " <<  _kbp.size2() << "\n" << flush; 
    _psi_AxB_dimer_basis.resize(0,0);
    _combAB.resize(0,0);
    _combA.resize(0,0);
    _combB.resize(0,0);
    // now the different spin types
            if (_doSinglets) {
                XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "   Evaluating singlets" << flush;
                // get singlet BSE Hamiltonian from _orbitalsAB
                Eigen::MatrixXd _Hamiltonian_AB = _eh_s.cast<double>();
                XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "   Setup Hamiltonian" << flush;
                const Eigen::MatrixXd _bseA_T = _orbitalsA->BSESingletCoefficients().block(0,0,_orbitalsA->BSESingletCoefficients().rows(),_levA).transpose().cast<double>();
                const Eigen::MatrixXd _bseB_T =_orbitalsB->BSESingletCoefficients().block(0,0,_orbitalsB->BSESingletCoefficients().rows(),_levB).transpose().cast<double>();
                
                JAB_singlet = ProjectExcitons(_bseA_T, _bseB_T, _Hamiltonian_AB);
                XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "   calculated singlet couplings " << flush;
            }



            if (_doTriplets) {
                XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "   Evaluating triplets" << flush;
                // get triplet BSE Hamiltonian from _orbitalsAB
                Eigen::MatrixXd _Hamiltonian_AB = _eh_s.cast<double>();
                XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "  Converted Hamiltonian to double" << flush;         
                const Eigen::MatrixXd _bseA_T = _orbitalsA->BSETripletCoefficients().block(0,0,_orbitalsA->BSETripletCoefficients().rows(),_levA).transpose().cast<double>();
                const Eigen::MatrixXd _bseB_T =_orbitalsB->BSETripletCoefficients().block(0,0,_orbitalsB->BSETripletCoefficients().rows(),_levB).transpose().cast<double>();
                JAB_triplet = ProjectExcitons(_bseA_T, _bseB_T, _Hamiltonian_AB);
                XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "   calculated triplet couplings " << flush;
            }
    
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "  Done with exciton couplings" << flush;
    return true;   
};


std::vector< Eigen::MatrixXd > BSECoupling::ProjectExcitons(const Eigen::MatrixXd& _bseA_T, const Eigen::MatrixXd& _bseB_T, 
                                  Eigen::MatrixXd& _H){
   
     // get projection of monomer excitons on dimer product functions
     Eigen::MatrixXd _proj_excA = _bseA_T* _kap;
     Eigen::MatrixXd _proj_excB =  _bseB_T* _kbp;
     
     _bseA_exc = _proj_excA.rows();
     _bseB_exc = _proj_excB.rows();
     _bse_exc=_bseA_exc+_bseB_exc;
   
     unsigned _ctAB=ctAB.rows();
     
     unsigned _ctBA=ctBA.rows();
     _ct=_ctAB+_ctBA;
     unsigned nobasisfunc=_H.rows();
     
     
     
     Eigen::MatrixXd fe_states=Eigen::MatrixXd::Zero(_bse_exc,nobasisfunc);
     fe_states.block(0,0, _bseA_exc, nobasisfunc )=_proj_excA;
     fe_states.block( _bseA_exc,0,_bseB_exc,nobasisfunc )=_proj_excB;
      
     Eigen::MatrixXd ct_states=Eigen::MatrixXd::Zero(_ct,nobasisfunc);
     
     //cout<< _ct<< "ct states"<<endl;
    if(_ct>0){ 
     //orthogonalize ct-states with respect to the FE states. 
       XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << " Orthogonalizing CT-states with respect to FE-states" << flush;
   
     if(_ctAB>0){
     ct_states.block( 0 ,0, _ctAB ,nobasisfunc)=ctAB;
    }
    if(_ctBA>0){
     ct_states.block( _ctAB, 0,_ctBA, nobasisfunc)=ctBA;
     }
       
       
        //orthogonalize ct-states with respect to FE states
     Eigen::MatrixXd correction=ct_states*fe_states.transpose()*fe_states;
      

    ct_states-=correction;    

    
     correction.resize(0,0);
     //normalize
    
     for (unsigned i=0;i<_ct;i++){
         double norm=ct_states.row(i).norm();
         if(norm<0.95){
            XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << " WARNING: CT-state "<< i<< " norm is only"<< norm << flush; 
         }
         ct_states.row(i)/=norm;
     }  
     }
      
     Eigen::MatrixXd projection(_bse_exc+_ct,nobasisfunc);
     
     
     XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << " merging projections into one vector  " << flush;
    
  projection.block(0 ,0, _bse_exc,nobasisfunc)=fe_states;
   
     if(_ct>0){
    projection.block(_bse_exc,0,_ct ,nobasisfunc )=ct_states;
     }
      XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "   Setting up coupling matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;
     // matrix _J
     
    //  E_A         J_AB        J_A_ABCT        J_A_BACT
    //  J_BA        E_B         J_B_ABCT        J_B_BACT
    //  J_ABCT_A    J_ABCT_B    E_ABCT          J_ABCT_BACT
    //  J_BACT_A   J_BACT_B    J_BACT_ABCT     E_BACT
     
     // I think this only works for hermitian/symmetric H so only in TDA
     // setup J
    
     Eigen::MatrixXd _J_dimer=projection*_H*projection.transpose();
      _H.resize(0,0);
   

    
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "   Setting up overlap matrix size "<< _bse_exc +_ct<<"x"<<_bse_exc +_ct << flush;
     // setup S
    
    Eigen::MatrixXd _S_dimer=projection*projection.transpose();
    
    projection.resize(0,0);
    if(tools::globals::verbose &&  _bse_exc+_ct<100){
         XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     XTP_LOG(xtp::logDEBUG, *_pLog) << "_J_dimer[Ryd]"<<flush;
     
     XTP_LOG(xtp::logDEBUG, *_pLog) << _J_dimer<<flush;
     XTP_LOG(xtp::logDEBUG, *_pLog) << "_S_dimer"<<flush;
     
     XTP_LOG(xtp::logDEBUG, *_pLog) << _S_dimer<<flush;
      XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    }
   
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_S_dimer);
    Eigen::MatrixXd Sm1=es.operatorInverseSqrt();
    _J_dimer=Sm1*_J_dimer*Sm1;
   
    
    if(tools::globals::verbose && _bse_exc+_ct<100){
         XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    XTP_LOG(xtp::logDEBUG, *_pLog) << "_J_ortho[Ryd]"<<flush;
    XTP_LOG(xtp::logDEBUG, *_pLog) << _J_dimer<<flush;
    XTP_LOG(xtp::logDEBUG, *_pLog) << "_S-1/2"<<flush;
    XTP_LOG(xtp::logDEBUG, *_pLog) << Sm1<<flush;
     XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
    }
     XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "   Smallest value of dimer overlapmatrix is "<<es.eigenvalues()(0)<< flush;
     
    std::vector< Eigen::MatrixXd >_J;
     
     XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "   Running Perturbation algorithm"<< flush;
    _J.push_back( Perturbation(_J_dimer));
    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp()  << "    Running Projection algorithm"<< flush;
    _J.push_back( Fulldiag(_J_dimer));
    
    
       if(tools::globals::verbose){
     XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     XTP_LOG(xtp::logDEBUG, *_pLog) << "Jeff_pert[Hrt]"<<flush;
     XTP_LOG(xtp::logDEBUG, *_pLog) << _J[0]<<flush;
     XTP_LOG(xtp::logDEBUG, *_pLog) << "Jeff_diag[Hrt]"<<flush;
     XTP_LOG(xtp::logDEBUG, *_pLog) << _J[1]<<flush;
     XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------"<<flush;
     }
      
     return _J;
}

Eigen::MatrixXd BSECoupling::Perturbation(const Eigen::MatrixXd& _J_dimer){
    
    Eigen::MatrixXd _J =Eigen::MatrixXd::Zero(_bse_exc, _bse_exc);
    bool _diag_ct = true;
    Eigen::MatrixXd _J_result=_J_dimer;
    if (_ct > 0 && _diag_ct) {

        Eigen::MatrixXd transformation = Eigen::MatrixXd::Identity(_bse_exc + _ct, _bse_exc + _ct);
        
        Eigen::MatrixXd Ct = _J_dimer.block(_bse_exc,_bse_exc,_ct, _ct);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ct);

        
        transformation.block(_bse_exc, _bse_exc ,_ct,_ct) = es.eigenvectors();

        Ct.resize(0, 0);

        if (tools::globals::verbose) {

            XTP_LOG(xtp::logDEBUG, *_pLog) << "FE state hamiltonian" << flush;
            XTP_LOG(xtp::logDEBUG, *_pLog) << _J_dimer.block(0,0, _bse_exc,_bse_exc) << flush;
            if (_ct > 0) {
                XTP_LOG(xtp::logDEBUG, *_pLog) << "eigenvalues of CT states" << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << es.eigenvalues() << flush;
            }

        }

     _J_result = transformation.transpose()*_J_dimer*transformation;
        if (tools::globals::verbose && _bse_exc + _ct < 100) {
            XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
            XTP_LOG(xtp::logDEBUG, *_pLog) << "_J_ortho[Hrt] CT-state diag" << flush;
            XTP_LOG(xtp::logDEBUG, *_pLog) << _J_result << flush;
            XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
        }
    }
    for (int stateA = 0; stateA < _levA; stateA++) {
        double Ea = _J_result(stateA, stateA);
        for (int stateB = 0; stateB < _levB; stateB++) {
            int stateBd = stateB + _bseA_exc;
            XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "   Calculating coupling between exciton A" << stateA + 1 << " and exciton B" << stateB + 1 << flush;
            double J = _J_result(stateA, stateBd);

            double Eb = _J_result(stateBd, stateBd);
            for (unsigned k = _bse_exc; k < (_bse_exc + _ct); k++) {
                double Eab = _J_result(k, k);
                if (std::abs(Eab - Ea) < 0.001) {
                    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "Energydifference between state A " << stateA + 1 << "and CT state " << k + 1 << " is " << Eab - Ea << "[Hrt]" << flush;
                }
                if (std::abs(Eab - Eb) < 0.001) {
                    XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "Energydifference between state B " << stateB + 1 << "and CT state " << k + 1 << " is " << Eab - Eb << "[Hrt]" << flush;

                }
                J += 0.5 * _J_result(k, stateA) * _J_result(k, stateBd)*(1 / (Ea - Eab) + 1 / (Eb - Eab)); // Have no clue why 0.5
            }
            _J(stateA, stateBd) = J;
            _J(stateBd, stateA) = J;


        }
    }

            
    return _J;
}


Eigen::MatrixXd BSECoupling::Fulldiag(const Eigen::MatrixXd& _J_dimer){
    Eigen::MatrixXd _J = Eigen::MatrixXd::Zero(_bse_exc, _bse_exc);
   
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_J_dimer);
  
    if (tools::globals::verbose && _bse_exc + _ct < 10) {
        XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
        XTP_LOG(xtp::logDEBUG, *_pLog) << "Eigenvectors of J" << flush;

        XTP_LOG(xtp::logDEBUG, *_pLog) << es.eigenvectors() << flush;
        XTP_LOG(xtp::logDEBUG, *_pLog) << "J_eigenvalues[Hrt]" << flush;
        XTP_LOG(xtp::logDEBUG, *_pLog) << es.eigenvalues() << flush;
        XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
    }
    //Calculate projection on subspace for every pair of excitons separately
    for (int stateA = 0; stateA < _levA; stateA++) {
        for (int stateB = 0; stateB < _levB; stateB++) {
            int stateBd = stateB + _bseA_exc;
            XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "   Calculating coupling between exciton A" << stateA + 1 << " and exciton B" << stateB + 1 << flush;
            std::vector<unsigned> index;
            std::vector<int> signvec;
            for (unsigned i = 0; i < _bse_exc + _ct; i++) {
                if (i == unsigned(stateA) || i == unsigned(stateBd)) {

                    double close = 0.0;
                    unsigned ind = 0;
                    int sign = 0;
                    //row
                    for (unsigned j = 0; j < _bse_exc + _ct; j++) {
                        bool check = true;
                        // if index i is already in index
                        // should not happen but if one vector was similar to two others.
                        for (unsigned l = 0; l < index.size(); l++) {
                            if (j == index[l]) {
                                check = false;
                                break;
                            }
                        }

                        if (check && std::abs(es.eigenvalues()(i, j)) > close) {
                            ind = j;
                            close = std::abs(es.eigenvalues()(i, j));
                            if (es.eigenvalues()(i, j) >= 0) {
                                sign = 1;
                            } else {
                                sign = -1;
                            }
                        }
                    }
                    index.push_back(ind);
                    signvec.push_back(sign);
                }
            }

            XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "   Order is: [Initial state n->nth eigenvalue]" << flush;
            XTP_LOG(xtp::logDEBUG, *_pLog) << "    A" << stateA + 1 << ":" << stateA + 1 << "->" << index[0] + 1 << " ";
            XTP_LOG(xtp::logDEBUG, *_pLog) << "    B" << stateB + 1 << ":" << stateBd + 1 << "->" << index[1] + 1 << " " << flush;

            //setting up transformation matrix _T and diagonal matrix _E for the eigenvalues;

            Eigen::MatrixXd _E = Eigen::MatrixXd::Zero(2, 2);
            Eigen::MatrixXd _T = Eigen::MatrixXd::Zero(2, 2);
            //find the eigenvectors which are most similar to the initial states

            //row 
            for (unsigned i = 0; i < 2; i++) {
                unsigned k = index[i];
                double sign = signvec[i];
                double normr = 1 / std::sqrt(es.eigenvectors()(stateA, k) * es.eigenvectors()(stateA, k) + es.eigenvectors()(stateBd, k) * es.eigenvectors()(stateBd, k));
                _T(0, i) = sign * es.eigenvectors()(stateA, k) * normr;
                _T(1, i) = sign * es.eigenvectors()(stateBd, k) * normr;
                _E(i, i) = es.eigenvectors()(k);
            }


            if ((_T(1, 1) * _T(0, 0) - _T(1, 0) * _T(0, 1)) < 0) {
                XTP_LOG(xtp::logDEBUG, *_pLog) << " Reduced state matrix is not in a right handed basis, multiplying second eigenvector by -1 " << flush;
                _T(0, 1) = -_T(0, 1);
                _T(1, 1) = -_T(1, 1);
            }

            if (tools::globals::verbose) {
                XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << "_T" << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << _T << flush;

            }

            Eigen::MatrixXd S_small = _T*_T.transpose();
            if (tools::globals::verbose) {

                XTP_LOG(xtp::logDEBUG, *_pLog) << "S_small" << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << S_small << flush;

            }
            //orthogonalize that matrix
            
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ss(S_small);
            Eigen::MatrixXd sm1=ss.operatorInverseSqrt();
            _E=sm1*_E*sm1;
            
            XTP_LOG(xtp::logDEBUG, *_pLog) << xtp::TimeStamp() << "   Smallest value of dimer overlapmatrix is " << ss.eigenvalues()(0) << flush;
            if (tools::globals::verbose) {

                XTP_LOG(xtp::logDEBUG, *_pLog) << "S-1/2" << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << sm1 << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << "E_ortho" << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << _E << flush;
            }
            _T = _T*sm1;
           
            if (tools::globals::verbose) {

                XTP_LOG(xtp::logDEBUG, *_pLog) << "T_ortho" << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << _T << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << "---------------------------------------" << flush;
            }


            Eigen::MatrixXd _J_small = _T*_E*_T.transpose();
            if (tools::globals::verbose) {
                XTP_LOG(xtp::logDEBUG, *_pLog) << "T_ortho*E_ortho*T_ortho^T" << flush;
                XTP_LOG(xtp::logDEBUG, *_pLog) << _J_small << flush;
            }

            _J(stateA, stateBd) = _J_small(0, 1);
            _J(stateBd, stateA) = _J_small(1, 0);

        }
    }
       
    return _J;
}


    
}}
