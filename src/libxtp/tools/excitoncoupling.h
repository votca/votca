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
    
class ExcitonCoupling : public  QMTool
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
    Logger      _log;

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
   
    _log.setReportLevel(  logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface( logINFO,    "\n... ...");
    _log.setPreface( logERROR,   "\n... ...");
    _log.setPreface( logWARNING, "\n... ...");
    _log.setPreface( logDEBUG,   "\n... ..."); 
    Property summary;
    Property& job_output = summary.add("output","");
    // get the corresponding object from the QMPackageFactory
    if(!_classical){
    Orbitals orbitalsA, orbitalsB, orbitalsAB;
    // load the QM data from serialized orbitals objects

    XTP_LOG( logDEBUG, _log) << " Loading QM data for molecule A from " << _orbA << flush;
    orbitalsA.ReadFromCpt(_orbA);
    
    XTP_LOG( logDEBUG, _log) << " Loading QM data for molecule B from " << _orbB << flush;
    orbitalsB.ReadFromCpt(_orbB);

    XTP_LOG( logDEBUG, _log) << " Loading QM data for dimer AB from " << _orbAB << flush;
    orbitalsAB.ReadFromCpt(_orbAB);
   
     BSECoupling bsecoupling; 
     bsecoupling.setLogger(&_log);
     bsecoupling.Initialize(_coupling_options);
  
     bsecoupling.CalculateCouplings( orbitalsA,orbitalsB, orbitalsAB );   
     std::cout << _log;

    Property& pair_summary = job_output.add("pair","");
    Property& type_summary = pair_summary.add("type","");
    bsecoupling.Addoutput(type_summary,orbitalsA,  orbitalsB);

    }
    
    else if (_classical){
        XTP_LOG( logDEBUG, _log) << "Calculating electronic coupling using classical transition charges." << _orbB << flush;
        std::vector< APolarSite*> seg1= APS_FROM_MPS(_mpsA, 0);
        std::vector< APolarSite*> seg2= APS_FROM_MPS(_mpsB, 0);
        
        PolarSeg Seg1 = PolarSeg(1,seg1);
        PolarSeg Seg2 = PolarSeg(2,seg2);
        XInteractor actor;
        actor.ResetEnergy();
        vec s = vec(0,0,0);
        double E = 0.0;
        for (APolarSite* site1:Seg1) {
          for (APolarSite* site2:Seg2) {            
            actor.BiasIndu(*site1, *site2, s);
            site1->Depolarize();
            site2->Depolarize();
            E += actor.E_f(*site1, *site2);             
          }
        }

    double J=E*conv::int2eV;  

    Property &pair_summary = job_output.add("pair","");
    pair_summary.setAttribute("idA", 1);
    pair_summary.setAttribute("idB", 2);
    pair_summary.setAttribute("typeA", _mpsA);
    pair_summary.setAttribute("typeB", _mpsB);
    Property & coupling_summary =pair_summary.add("Coupling",""); 
    coupling_summary.setAttribute("jABstatic", J);
    }
    
    tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 1, "");
     
    std::ofstream ofs (_output_file.c_str(), std::ofstream::out);
    ofs << job_output;    
    ofs.close();
    return true;
}


}}


#endif
