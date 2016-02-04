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

#ifndef _VOTCA_XTP_EXCITONCOUPLINGH_H
#define _VOTCA_XTP_EXCITONCOUPLINGH_H

#include <stdio.h>

#include <votca/xtp/logger.h>
#include <votca/tools/constants.h>
#include <votca/xtp/bsecoupling.h>
#include <votca/xtp/qmpackagefactory.h>

namespace votca { namespace xtp {
    using namespace std;
    
class ExcitonCoupling : public QMTool
{
public:

    ExcitonCoupling() { };
   ~ExcitonCoupling() { };

    string Identify() { return "excitoncoupling"; }

    void   Initialize(Property *options);
    bool   Evaluate();

 

private:
    
    string      _orbA, _orbB, _orbAB;
    string      _logA, _logB, _logAB;
    int         _levA, _levB;
   // int         _trimA, _trimB;
    double      _degeneracy;

    string      _spintype;
    Property    _package_options; 
    
    string      _output_file;
    bool        _classical;
    string      _mpsA;
    string      _mpsB;
    
    
    Logger      _log;

};

void ExcitonCoupling::Initialize(Property* options) 
{

   // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "xtp" );
    std::string key = "options." + Identify();  
    _classical=false;
    if ( options->exists(key+".classical")) {
        _classical = options->get(key+".classical").as<bool>();
        }
    
    if(!_classical){
    
        _degeneracy = options->get(key + ".degeneracy").as<double> ();
        _spintype   = options->get(key + ".type").as<string> ();



        _orbA  = options->get(key + ".moleculeA.orbitals").as<string> ();
        _orbB  = options->get(key + ".moleculeB.orbitals").as<string> ();
        _orbAB = options->get(key + ".dimerAB.orbitals").as<string> ();

        _levA  = options->get(key + ".moleculeA.states").as<int> ();
        _levB  = options->get(key + ".moleculeB.states").as<int> ();
    }
    else{
        _mpsA= options->get(key + ".moleculeA.mps").as<string> ();
        _mpsB= options->get(key + ".moleculeB.mps").as<string> ();
    }
    _output_file = options->get(key + ".output").as<string> ();

    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    
}

bool ExcitonCoupling::Evaluate() {
    Property *_job_output=NULL;
    Property _summary; 
    _log.setReportLevel( logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(logINFO,    "\n... ...");
    _log.setPreface(logERROR,   "\n... ...");
    _log.setPreface(logWARNING, "\n... ...");
    _log.setPreface(logDEBUG,   "\n... ..."); 

    // get the corresponding object from the QMPackageFactory
    if(!_classical){
    Orbitals _orbitalsA, _orbitalsB, _orbitalsAB;
    // load the QM data from serialized orbitals objects

    std::ifstream ifa( (_orbA ).c_str());
    LOG(logDEBUG, _log) << " Loading QM data for molecule A from " << _orbA << flush;
    boost::archive::binary_iarchive ia(ifa);
    ia >> _orbitalsA;
    ifa.close();
    
    std::ifstream ifb( (_orbB ).c_str());
    LOG(logDEBUG, _log) << " Loading QM data for molecule B from " << _orbB << flush;
    boost::archive::binary_iarchive ib(ifb);
    ib >> _orbitalsB;
    ifb.close();
    
    
    std::ifstream ifab( (_orbAB ).c_str());
    LOG(logDEBUG, _log) << " Loading QM data for dimer AB from " << _orbAB << flush;
    boost::archive::binary_iarchive iab(ifab);
    iab >> _orbitalsAB;
    ifab.close();

     BSECoupling _bsecoupling; 
     _bsecoupling.setLogger(&_log);
          
     ub::matrix<float> _JAB_singlet;
     ub::matrix<float> _JAB_triplet;

     //bool _calculate_integrals = _bsecoupling.CalculateCouplings_OLD( &_orbitalsA, &_orbitalsB, &_orbitalsAB, &_JAB_singlet );   
     bool _calculate_integrals = _bsecoupling.CalculateCouplings( &_orbitalsA, &_orbitalsB, &_orbitalsAB, &_JAB_singlet, &_JAB_triplet, _spintype );   
     std::cout << _log;
 
    if ( _calculate_integrals ){ 
     // output the results
    
    _job_output = &_summary.add("output","");
    Property *_pair_summary = &_job_output->add("pair","");
    Property *_type_summary = &_pair_summary->add("type","");
    if ( _spintype == "singlets" || _spintype == "all" ){
        Property *_singlet_summary = &_type_summary->add("singlets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB < _levB  ; ++stateB ) {
               float JAB = _bsecoupling.getSingletCouplingElement( stateA , stateB, &_orbitalsA, &_orbitalsB, &_JAB_singlet, _degeneracy );
               Property *_coupling_summary = &_singlet_summary->add("coupling", boost::lexical_cast<string>(JAB)); 
               float energyA = _orbitalsA.BSESingletEnergies()[stateA]*27.21138386/2.0;
               float energyB = _orbitalsB.BSESingletEnergies()[stateB]*27.21138386/2.0;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", energyA);
               _coupling_summary->setAttribute("energyB", energyB);
           } 
        }
    }
    if ( _spintype == "triplets" || _spintype == "all" ){
        Property *_triplet_summary = &_type_summary->add("triplets","");
        for (int stateA = 0; stateA < _levA ; ++stateA ) {
           for (int stateB = 0; stateB < _levB  ; ++stateB ) {
               float JAB = _bsecoupling.getTripletCouplingElement( stateA , stateB, &_orbitalsA, &_orbitalsB, &_JAB_triplet, _degeneracy );
               Property *_coupling_summary = &_triplet_summary->add("coupling", boost::lexical_cast<string>(JAB)); 
               float energyA = _orbitalsA.BSETripletEnergies()[stateA]*27.21138386/2.0;
               float energyB = _orbitalsB.BSETripletEnergies()[stateB]*27.21138386/2.0;
               _coupling_summary->setAttribute("excitonA", stateA);
               _coupling_summary->setAttribute("excitonB", stateB);
               _coupling_summary->setAttribute("energyA", energyA);
               _coupling_summary->setAttribute("energyB", energyB);
           } 
        }
    }    
    
   
   
    }
    }
    
    else if (_classical){
        std::vector<APolarSite*> seg1=APS_FROM_MPS(_mpsA, 0);
        std::vector<APolarSite*> seg2=APS_FROM_MPS(_mpsB, 0);
        
        PolarSeg* Seg1 = new PolarSeg(1,seg1);
        PolarSeg* Seg2 = new PolarSeg(2,seg2);
        XInteractor actor;
        actor.ResetEnergy();
        Seg1->CalcPos();
        Seg2->CalcPos();
        vec s = Seg1->getPos() - Seg2->getPos();
        //LOG(logINFO, *pLog) << "Evaluate pair for debugging " << Seg1->getId() << ":" <<Seg2->getId() << " Distance "<< abs(s) << flush; 
        PolarSeg::iterator pit1;
        PolarSeg::iterator pit2;
        double E = 0.0;
        for (pit1 = Seg1->begin(); pit1 < Seg1->end(); ++pit1) {
            for (pit2 = Seg2->begin(); pit2 < Seg2->end(); ++pit2) {
                actor.BiasIndu(*(*pit1), *(*pit2), s);
                (*pit1)->Depolarize();
                (*pit2)->Depolarize();
                E += actor.E_f(*(*pit1), *(*pit2));
                               
                }
            }
    double J=E*conv::int2eV;  
    
    _job_output = &_summary.add("output","");
    Property *_pair_summary = &_job_output->add("pair","");
    _pair_summary->setAttribute("idA", 1);
    _pair_summary->setAttribute("idB", 2);
    _pair_summary->setAttribute("typeA", _mpsA);
    _pair_summary->setAttribute("typeB", _mpsB);
    Property *_coupling_summary = &_pair_summary->add("Coupling",""); 
    _coupling_summary->setAttribute("jABstatic", J);
    }
    
    votca::tools::PropertyIOManipulator iomXML(votca::tools::PropertyIOManipulator::XML, 1, "");
     
    std::ofstream ofs (_output_file.c_str(), std::ofstream::out);
    ofs << *_job_output;    
    ofs.close();
    return true;
}


}}


#endif