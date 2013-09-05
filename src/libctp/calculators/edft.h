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


#ifndef _CALC_DFT_ENERGIES_H
#define	_CALC_DFT_ENERGIES_H

#include <votca/ctp/segment.h>
#include <votca/ctp/orbitals.h>

#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/parallelxjobcalc.h>
#include <unistd.h>

#include <fstream>
#include <sys/stat.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>

namespace votca { namespace ctp {

/**
* \brief Site energies and orbitals for QM molecules
*
* QM orbitals and energies for all molecules
* Requires a first-principles package, i.e. GAUSSIAN installation
*
* Callname: edft
*/

class EDFT : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{
public:

    EDFT() {};
   ~EDFT() {};

    string  Identify() { return "edft"; }
    void    Initialize(Topology *top, Property *options);
    void    ParseOrbitalsXML(Topology *top, Property *options);
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *thread);

    void    CleanUp();


private:

    // what to do
    bool                _do_input;
    bool                _do_run;
    bool                _do_parse;

    // what to write in the storage
    bool                _store_orbitals;

    //bool   _maverick;
    string _outParent;
    string _outMonDir;
    
    string _package;
    Property _package_options;   

};

void EDFT::CleanUp() {

}

void EDFT::Initialize(Topology *top, Property *options) {

    _do_input = false;
    _do_run = false;
    _do_parse = false;

    _maverick = (_nThreads == 1) ? true : false;
    
    string key = "options." + Identify();
    _jobfile = options->get(key + ".job_file").as<string>();

    string _package_xml = options->get(key+".package").as<string> ();

    string _tasks_string = options->get(key+".tasks").as<string> ();
    if (_tasks_string.find("input") != std::string::npos) _do_input = true;
    if (_tasks_string.find("run") != std::string::npos) _do_run = true;
    if (_tasks_string.find("parse") != std::string::npos) _do_parse = true;    

    string _store_string = options->get(key+".store").as<string> ();
    if (_store_string.find("orbitals") != std::string::npos) _store_orbitals = true;
    
    load_property_from_xml( _package_options, _package_xml.c_str() );    
    key = "package";
    _package = _package_options.get(key+".name").as<string> ();
    
    // register all QM packages (Gaussian, turbomole, nwchem))
    QMPackageFactory::RegisterAll(); 

}


Job::JobResult EDFT::EvalJob(Topology *top, Job *job, QMThread *opThread) {

    string output;
    
    bool _run_status;
    bool _parse_log_status;
    bool _parse_orbitals_status;

    FILE *out;
    Orbitals _orbitals;
    Job::JobResult jres = Job::JobResult();

    vector < Segment* > segments;    
    Segment *seg = top->getSegment( boost::lexical_cast<int>( job->getTag() ));
    segments.push_back( seg );

    Logger* pLog = opThread->getLogger();
    LOG(logINFO,*pLog) << TimeStamp() << " Evaluating site " << seg->getId() << flush; 

    // log, com, and orbital files will be stored in ORB_FILES/package_name/frame_x/mol_ID/
    // extracted information will be stored in  ORB_FILES/molecules/frame_x/molecule_ID.orb
    
    string edft_work_dir = "OR_FILES";

    string orbitals_storage_dir = "molecules";
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());      
    
    string ID   = boost::lexical_cast<string>( seg->getId() );
    string qmpackage_work_dir  = edft_work_dir + "/" + _package + "/" + frame_dir + "/mol_" + ID;
    string ORB_DIR = edft_work_dir + "/" + orbitals_storage_dir + "/" + frame_dir;

    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
     
    boost::filesystem::create_directories( qmpackage_work_dir );     
    boost::filesystem::create_directories(ORB_DIR);        
    
   _qmpackage->setLog( pLog );  
   
   _qmpackage->Initialize( &_package_options );

   _qmpackage->setRunDir( qmpackage_work_dir );
  

    // if asked, prepare the input files
    if ( _do_input ) {
        _qmpackage->WriteInputFile( segments );
    }
        
   // Run the executable
    if ( _do_run ) {
        _run_status = _qmpackage->Run( );
        if ( !_run_status ) {
            output += "run failed; " ;
            LOG(logERROR,*pLog) << _package << " run failed" << flush;
            jres.setOutput( output ); 
            jres.setStatus(Job::FAILED);
            delete _qmpackage;
            return jres;
        } else {
            output += "run completed; " ;
        }
    }
        
   // Parse log files
    if ( _do_parse ) {
        _parse_log_status = _qmpackage->ParseLogFile( &_orbitals );
        if ( !_parse_log_status ) {
            output += "log incomplete; ";
            LOG(logERROR,*pLog) << "QM log incomplete" << flush;
            jres.setOutput( output ); 
            jres.setStatus(Job::FAILED);
            delete _qmpackage;
            return jres;
        } else {
            output += "log parsed; " ;
        }

       // Parse orbitals file
       _parse_orbitals_status = _qmpackage->ParseOrbitalsFile( &_orbitals );
        if ( !_parse_orbitals_status ) {
            output += "orbitals failed; " ;
            LOG(logERROR,*pLog) << "QM orbitals not parsed" << flush;
            jres.setOutput( output ); 
            jres.setStatus(Job::FAILED);
            delete _qmpackage;
            return jres;
        } else {
            output += "orbitals parsed; " ;
        }
    }
   
   // Clean run
   _qmpackage->CleanUp();
   delete _qmpackage;
        
    // GENERATE OUTPUT AND FORWARD TO PROGRESS OBSERVER (RETURN)
    jres.setStatus(Job::COMPLETE);


    // save orbitals
    string ORB_FILE = "molecule_" + ID + ".orb";
    LOG(logDEBUG,*pLog) << "Serializing to " <<  ORB_FILE << flush;
    std::ofstream ofs( (ORB_DIR + "/" + ORB_FILE).c_str() );
    boost::archive::binary_oarchive oa( ofs );
    oa << _orbitals;
    ofs.close();
                                
    LOG(logINFO,*pLog) << TimeStamp() << " Finished evaluating site " << seg->getId() << flush; 
    
    // output of the JOB 
    jres.setOutput( output );

    // dump the LOG
    cout << *pLog;
    
    return jres;
    
   
}

}}

#endif	/* _CALC_DFT_ENERGIES_H */
