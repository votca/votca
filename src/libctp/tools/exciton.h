/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef _VOTCA_CTP_EXCITON_H
#define _VOTCA_CTP_EXCITON_H

#include <stdio.h>

#include <votca/ctp/logger.h>
#include <votca/ctp/gwbse.h>
#include <votca/ctp/qmpackagefactory.h>

namespace votca { namespace ctp {
    using namespace std;
    
class Exciton : public QMTool
{
public:

    Exciton() { };
   ~Exciton() { };

    string Identify() { return "exciton"; }

    void   Initialize(Property *options);
    bool   Evaluate();

    GWBSE _gwbse;
 

private:
    
    string      _orbfile;
    string      _logfile;

    string      _package;
    Property    _package_options;
    Property    _gwbse_options;
    
    string      _output_file;
    
    Logger      _log;

};

void Exciton::Initialize(Property* options) {

            // update options with the VOTCASHARE defaults   
            UpdateWithDefaults( options );

            string key = "options." + Identify();
           _output_file =  options->get(key + ".archive").as<string>();

           // options for GWBSE package
           string _gwbse_xml = options->get(key+".gwbse").as<string> ();
           //cout << endl << "... ... Parsing " << _package_xml << endl ;
           load_property_from_xml( _gwbse_options, _gwbse_xml.c_str() );    
           
           
            // something unique
           _package = options->get(key + ".package").as<string> ();
       
           // orbitals file or pure DFT output
           _orbfile  = options->get(key + ".molecule.orbitals").as<string> ();
           _logfile  = options->get(key + ".molecule.log").as<string> ();

            
    
  
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/qmpackages/") + _package + string("_idft_pair.xml");
    load_property_from_xml( _package_options, xmlFile );    

    // register all QM packages (Gaussian, TURBOMOLE, etc)
    QMPackageFactory::RegisterAll();
    
}

bool Exciton::Evaluate() {

    _log.setReportLevel( logDEBUG );
    _log.setMultithreading( true );
    
    _log.setPreface(logINFO,    "\n... ...");
    _log.setPreface(logERROR,   "\n... ...");
    _log.setPreface(logWARNING, "\n... ...");
    _log.setPreface(logDEBUG,   "\n... ..."); 

    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
   _qmpackage->setLog( &_log );       
   _qmpackage->Initialize( &_package_options );

     Orbitals _orbitals;
  
    _qmpackage->setOrbitalsFileName( _orbfile );
    int _parse_orbitals_status = _qmpackage->ParseOrbitalsFile( &_orbitals );

    _qmpackage->setLogFileName( _logfile );
    int _parse_log_status = _qmpackage->ParseLogFile( &_orbitals );
    
    _gwbse.setLogger(&_log);
    _gwbse.Initialize( &_gwbse_options );
     bool _evaluate = _gwbse.Evaluate( &_orbitals );
     std::cout << _log;
 
     
     // output the results
/*    Property _summary; 
    Property *_job_output = &_summary.add("output","");
    Property *_pair_summary = &_job_output->add("pair","");
    int HOMO_A = _orbitalsA.getNumberOfElectrons();
    int HOMO_B = _orbitalsB.getNumberOfElectrons();
    int LUMO_A = HOMO_A + 1;
    int LUMO_B = HOMO_B + 1;
    _pair_summary->setAttribute("homoA", HOMO_A);
    _pair_summary->setAttribute("homoB", HOMO_B);
    for (int levelA = HOMO_A - _levA +1; levelA <= LUMO_A + _levA - 1; ++levelA ) {
        for (int levelB = HOMO_B - _levB + 1; levelB <= LUMO_B + _levB -1 ; ++levelB ) {        
                double JAB = _overlap.getCouplingElement( levelA , levelB, &_orbitalsA, &_orbitalsB, &_JAB, _degeneracy );
                Property *_overlap_summary = &_pair_summary->add("overlap", boost::lexical_cast<string>(JAB)); 
                double energyA = _orbitalsA.getEnergy( levelA );
                double energyB = _orbitalsB.getEnergy( levelB );
                _overlap_summary->setAttribute("orbA", levelA);
                _overlap_summary->setAttribute("orbB", levelB);
                //_overlap_summary->setAttribute("jAB", JAB);
                _overlap_summary->setAttribute("eA", energyA);
                _overlap_summary->setAttribute("eB", energyB);
        }
    }
 */
     
           
            
       LOG(logDEBUG,_log) << "Saving data to " << _output_file << flush;
       std::ofstream ofs( ( _output_file).c_str() );
       boost::archive::binary_oarchive oa( ofs );


       
       oa << _orbitals;
       ofs.close();

     
     
     
     
    Property _summary; 
    Property *_job_output = &_summary.add("output","");
    votca::tools::PropertyIOManipulator iomXML(votca::tools::PropertyIOManipulator::XML, 1, "");
     
    //ofs (_output_file.c_str(), std::ofstream::out);
    //ofs << *_job_output;    
    //ofs.close();
    
    return true;
}


}}


#endif
