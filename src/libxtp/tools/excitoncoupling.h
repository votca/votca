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

#ifndef _VOTCA_XTP_EXCITONCOUPLINGH_H
#define _VOTCA_XTP_EXCITONCOUPLINGH_H


#include <votca/xtp/qmtool.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/xinteractor.h>

#include <stdio.h>
#include <votca/tools/constants.h>

#include <votca/tools/constants.h>
#include <votca/xtp/bsecoupling.h>

#include <votca/xtp/qmpackagefactory.h>

namespace votca { namespace xtp {
    using namespace std;
    
class ExcitonCoupling : public  xtp::QMTool
{
public:

    ExcitonCoupling() { };
   ~ExcitonCoupling() { };

    string Identify() { return "excitoncoupling"; }

    void   Initialize(Property *options);
    bool   Evaluate();

 

private:
    
    string      _orbA, _orbB, _orbAB;
   // int         _trimA, _trimB;
    
    Property    _coupling_options; 
    
    string      _output_file;
    bool        _classical;
    //bool        _doSinglets;
    //bool        _doTriplets;
    string      _mpsA;
    string      _mpsB;  
    xtp::Logger      _log;

};

void ExcitonCoupling::Initialize(Property* options) 
{
   // _doSinglets=false;
   // _doTriplets=false;
   // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options, "xtp" );
    std::string key = "options." + Identify();  
    _classical=false;
    if ( options->exists(key+".classical")) {
        _classical = options->get(key+".classical").as<bool>();
        
        }
    else{
        _classical=false;
    }
    
    if(!_classical){
        
        string _coupling_xml=options->get(key + ".bsecoupling_options").as<string>();
        load_property_from_xml(_coupling_options, _coupling_xml.c_str());
        
        _orbA  = options->get(key + ".orbitalsA").as<string> ();
        _orbB  = options->get(key + ".orbitalsB").as<string> ();
        _orbAB = options->get(key + ".orbitalsAB").as<string> ();

    }
    else{
        _mpsA= options->get(key + ".mpsA").as<string> ();
        _mpsB= options->get(key + ".mpsB").as<string> ();
    }
    _output_file = options->get(key + ".output").as<string> ();

    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    
}

bool ExcitonCoupling::Evaluate() {
    Property *_job_output=NULL;
    Property _summary; 
    _log.setReportLevel(  xtp::logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface( xtp::logINFO,    "\n... ...");
    _log.setPreface( xtp::logERROR,   "\n... ...");
    _log.setPreface( xtp::logWARNING, "\n... ...");
    _log.setPreface( xtp::logDEBUG,   "\n... ..."); 

    // get the corresponding object from the QMPackageFactory
    if(!_classical){
    Orbitals _orbitalsA, _orbitalsB, _orbitalsAB;
    // load the QM data from serialized orbitals objects

    XTP_LOG( xtp::logDEBUG, _log) << " Loading QM data for molecule A from " << _orbA << flush;
    _orbitalsA.ReadFromCpt(_orbA);
    
    XTP_LOG( xtp::logDEBUG, _log) << " Loading QM data for molecule B from " << _orbB << flush;
    _orbitalsB.ReadFromCpt(_orbB);

    XTP_LOG( xtp::logDEBUG, _log) << " Loading QM data for dimer AB from " << _orbAB << flush;
    _orbitalsAB.ReadFromCpt(_orbAB);
   
     BSECoupling _bsecoupling; 
     _bsecoupling.setLogger(&_log);
     _bsecoupling.Initialize(&_coupling_options);
     //bool _doSinglets=_bsecoupling.get_doSinglets();   
     //bool _dotriplets=_bsecoupling.get_doTriplets();   
     

     //bool _calculate_integrals = _bsecoupling.CalculateCouplings_OLD( &_orbitalsA, &_orbitalsB, &_orbitalsAB, &_JAB_singlet );   
     bool _calculate_integrals = _bsecoupling.CalculateCouplings( &_orbitalsA, &_orbitalsB, &_orbitalsAB );   
     std::cout << _log;
 
    if ( _calculate_integrals ){ 
     // output the results
    
    _job_output = &_summary.add("output","");
    Property *_pair_summary = &_job_output->add("pair","");
    Property *_type_summary = &_pair_summary->add("type","");
    _bsecoupling.addoutput(_type_summary,&_orbitalsA, & _orbitalsB);
    
   
   
    }
    }
    
    else if (_classical){
        XTP_LOG( xtp::logDEBUG, _log) << "Calculating electronic coupling using classical transition charges." << _orbB << flush;
        std::vector< xtp::APolarSite*> seg1= xtp::APS_FROM_MPS(_mpsA, 0);
        std::vector< xtp::APolarSite*> seg2= xtp::APS_FROM_MPS(_mpsB, 0);
        
        xtp::PolarSeg* Seg1 = new  xtp::PolarSeg(1,seg1);
         xtp::PolarSeg* Seg2 = new  xtp::PolarSeg(2,seg2);
         xtp::XInteractor actor;
        actor.ResetEnergy();
        vec s = vec(0,0,0);
        
        //XTP_LOG(logINFO, *pLog) << "Evaluate pair for debugging " << Seg1->getId() << ":" <<Seg2->getId() << " Distance "<< abs(s) << flush; 
         xtp::PolarSeg::iterator pit1;
         xtp::PolarSeg::iterator pit2;
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
    
    tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1, "");
     
    std::ofstream ofs (_output_file.c_str(), std::ofstream::out);
    ofs << *_job_output;    
    ofs.close();
    return true;
}


}}


#endif
