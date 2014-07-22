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


#ifndef _VOTCA_CTP_DMA_H
#define	_VOTCA_CTP_DMA_H

#include <votca/ctp/segment.h>

#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/parallelxjobcalc.h>
#include <unistd.h>

#include <fstream>
#include <sys/stat.h>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

using boost::format;

namespace votca { namespace ctp {

/**
* \brief Distributed multipole analysis using Gaussian input
*
* Evaluates distributed multipoles
* Requires GAUSSIAN and GDMA
*
* Callname: dma
*/

class DMA : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{
public:

    DMA() {};
   ~DMA() {};

    string   Identify() { return "dma"; }
    void     Initialize(Property *options);
    void     WriteJobFile(Topology *top);
    
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *thread);

private:

    // what to do
    bool                _do_input;
    bool                _do_orbitals;
    bool                _do_dma;

    string _jobFile;
            
    string _package;
    Property _package_options;   

};

void DMA::Initialize(Property *options) {

    _do_input = false;
    _do_orbitals = false;
    _do_dma = false;
           
    _maverick = (_nThreads == 1) ? true : false;
    
    string key = "options." + Identify();
    string _package_xml = options->get(key+".package").as<string> ();
    
    string _tasks_string = options->get(key+".tasks").as<string> ();
    if (_tasks_string.find("input") != std::string::npos) _do_input = true;
    if (_tasks_string.find("orbitals") != std::string::npos) _do_orbitals = true;
    if (_tasks_string.find("dma") != std::string::npos) _do_dma = true;
    
    key = "options." + Identify() +".job";
    _jobfile = options->get(key + ".file").as<string>();

    load_property_from_xml( _package_options, _package_xml.c_str() );    
    key = "package";
    _package = _package_options.get(key+".name").as<string> ();
    if ( _package != "gaussian" ) { throw runtime_error("\nERROR: Package: " + _package + "not supported. Use Gaussian for GDMA"); }
    
    // register all QM packages (Gaussian, Turbomole, NWChem))
    QMPackageFactory::RegisterAll(); 

}

void DMA::WriteJobFile(Topology *top) {

    cout << endl << "... ... Writing job file: " << flush;
    ofstream ofs;
    ofs.open(_jobfile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);
 
    ofs << "<jobs>" << endl;   

    std::vector<Segment*> &segments = top->Segments();
    std::vector<Segment*>::iterator sit;
 
    int jobCount = 0;

    for (sit = segments.begin(); sit != segments.end(); ++sit) {
    
        int id = ++jobCount;
        string tag = "";

        Property Input;
        Property *pInput = &Input.add("input","");
        Property *pSegment =  &pInput->add("segment" , "" ) ;
        //pSegment->setAttribute<string>("type", sit->getName() );
        //pSegment->setAttribute<int>("id", sit->getId() );
        Job job(id, tag, Input, Job::AVAILABLE );
        job.ToStream(ofs,"xml");
    }
     

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    
    cout << jobCount << " jobs" << flush;
    
}


Job::JobResult DMA::EvalJob(Topology *top, Job *job, QMThread *opThread) {

    string output;
    
    bool _orbitals_status;
    bool _formchk_status;
    bool _dma_status;


    Job::JobResult jres = Job::JobResult();
    Property _job_input = job->getInput();  
    list<Property*> lSegments = _job_input.Select( "segment" );  
    
    vector < Segment* > segments;    
    int segId = lSegments.front()->getAttribute<int>( "id" );
    string segType = lSegments.front()->getAttribute<string>( "type" );
    
    Segment *seg = top->getSegment( segId );
    assert( seg->getName() == segType ); 
    segments.push_back( seg );
    Logger* pLog = opThread->getLogger();
    LOG(logINFO,*pLog) << TimeStamp() << " Evaluating site " << seg->getId() << flush; 

    // log, com, and orbital files will be stored in ORB_FILES/package_name/frame_x/mol_ID/
    // extracted information will be stored in  ORB_FILES/molecules/frame_x/molecule_ID.orb
    
    string dma_work_dir = "DMA_FILES";
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());      
    string ID   = boost::lexical_cast<string>( seg->getId() );

    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
    
   _qmpackage->setLog( pLog );  
   _qmpackage->Initialize( &_package_options );

   string qmpackage_work_dir;
   
   qmpackage_work_dir  = dma_work_dir + "/" + _package + "/" + frame_dir + "/mol_" + ID;
   boost::filesystem::create_directories( qmpackage_work_dir );     
    
   _qmpackage->setRunDir( qmpackage_work_dir );
   
    // if asked, prepare the input files
    if ( _do_input ) {
        _qmpackage->WriteInputFile( segments );
    }
        
   // Run the Gaussian executable
    if ( _do_orbitals ) {
        _orbitals_status = _qmpackage->Run( );
        if ( !_orbitals_status ) {
            output += "run failed; " ;
            LOG(logERROR,*pLog) << _package << " run failed" << flush;
            jres.setOutput( output ); 
            jres.setStatus(Job::FAILED);
            delete _qmpackage;
            return jres;
        } else {
            output += "orbitals completed; " ;
        }
    }
        
   // Parse log files
    if ( _do_dma ) {
        _dma_status = true;
        if ( !_dma_status ) {
            output += "DMA job incomplete; ";
            LOG(logERROR,*pLog) << "DMA job incomplete" << flush;
            jres.setOutput( output ); 
            jres.setStatus(Job::FAILED);
           return jres;
        } else {
            output += "DMA complete; " ;
        }
    }

   
   // Clean run
   _qmpackage->CleanUp();
   delete _qmpackage;
        
    LOG(logINFO,*pLog) << TimeStamp() << " Finished evaluating site " << seg->getId() << flush; 
 
    Property _job_summary;
        Property *_output_summary = &_job_summary.add("output","");
        Property *_segment_summary = &_output_summary->add("segment","");
         string segName = seg->getName();
         segId = seg->getId();
        _segment_summary->setAttribute("id", segId);
    
    // output of the JOB 
    jres.setOutput( _job_summary );
    jres.setStatus(Job::COMPLETE);

    
    return jres;

}

}}

#endif	/* _VOTCA_CTP_DMA_H */
