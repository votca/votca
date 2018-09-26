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


#ifndef _VOTCA_XTP_EDFT_H
#define	_VOTCA_XTP_EDFT_H

#include <votca/xtp/segment.h>
#include <votca/xtp/orbitals.h>

#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/parallelxjobcalc.h>
#include <unistd.h>

#include <fstream>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

using boost::format;

namespace votca { namespace xtp {

/**
* \brief Site energies and orbitals for QM molecules
*
* Evaluates orbitals and energies for all molecules
* Requires a first-principles package, i.e. GAUSSIAN, TURBOMOLE, NWChem
*
* Callname: edft
*/

class EDFT : public ParallelXJobCalc< vector<Job*>,Job*, Job::JobResult >
{
public:

    EDFT() {};
   ~EDFT() {};

    string   Identify() { return "edft"; }
    void     Initialize(Property *options);
    void     WriteJobFile(Topology *top);
    
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *thread);

private:

    // what to do
    bool                _do_input;
    bool                _do_run;
    bool                _do_parse;
    bool                _do_trim;
 
    // what to write in the storage
    bool                _store_orbitals;


    string _outParent;
    string _jobFile;
            
    string _package;
    Property _package_options;   
    


};

void EDFT::Initialize(Property *options) {

    _do_input = false;
    _do_run = false;
    _do_parse = false;
    _do_trim = false;
    
 
        
    _maverick = (_nThreads == 1) ? true : false;
    
    string key = "options." + Identify();
    string _package_xml = options->get(key+".dftpackage").as<string> ();
    
    string _tasks_string = options->get(key+".tasks").as<string> ();
    if (_tasks_string.find("input") != std::string::npos) _do_input = true;
    if (_tasks_string.find("run") != std::string::npos) _do_run = true;
    if (_tasks_string.find("trim") != std::string::npos) _do_trim = true;
    if (_tasks_string.find("parse") != std::string::npos) _do_parse = true;    
   
    
    string _store_string = options->get(key+".store").as<string> ();
    if (_store_string.find("orbitals") != std::string::npos) _store_orbitals = true;
  
    
    key = "options."+Identify();

        if ( options->exists(key+".job_file")) {
            _jobfile = options->get(key+".job_file").as<string>();
        }
        else {
            throw std::runtime_error("Job-file not set. Abort.");
        }

    
    load_property_from_xml( _package_options, _package_xml.c_str() );    
    key = "package";
    _package = _package_options.get(key+".name").as<string> ();


  
    
    
    // register all QM packages (Gaussian, Turbomole, NWChem))
    QMPackageFactory::RegisterAll(); 

}

void EDFT::WriteJobFile(Topology *top) {

    cout << endl << "... ... Writing job file: " << flush;
    ofstream ofs;
    ofs.open(_jobfile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);
 
    ofs << "<jobs>" << endl;   

    QMNBList::iterator pit;
    QMNBList &nblist = top->NBList();    
    
            
    int jobCount = 0;
    if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
    } 

    // regenerate the list of bridging segments for every pair 
    // (Donor - Bridge1 - Bridge2 - ... - Acceptor) type
    nblist.GenerateSuperExchange();
    
    map< int,Segment* > segments;
    map< int,Segment* >::iterator sit;

    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {
        
        int id1 = (*pit)->Seg1()->getId();
        int id2 = (*pit)->Seg2()->getId();
	segments[id1] = (*pit)->Seg1();
        segments[id2] = (*pit)->Seg2();
        
        /* loop over bridging segments if any and add them to the map
           this in principle is not needed since all pairs between 
           donors, acceptors, and bridges are already in the list 
         */
        vector<Segment*> bridges = (*pit)->getBridgingSegments();
        for ( vector<Segment*>::const_iterator bsit = bridges.begin(); bsit != bridges.end(); bsit++ ) {
            //cout << "Bridging segment " << (*bsit)->getId() << " : " <<  (*bsit)->getName() << endl;
            segments[ (*bsit)->getId() ] = (*bsit);
        }

    }
    

    
    for (sit = segments.begin(); sit != segments.end(); ++sit) {
    
        int id = ++jobCount;
        string tag = "";

        Property Input;
        Property *pInput = &Input.add("input","");
        Property *pSegment =  &pInput->add("segment" , (format("%1$s") % sit->first).str() );
        pSegment->setAttribute<string>("type", sit->second->getName() );
        pSegment->setAttribute<int>("id", sit->second->getId() );
        Job job(id, tag, Input, Job::AVAILABLE );
        job.ToStream(ofs,"xml");
    }
     

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    
    cout << jobCount << " jobs" << flush;
    
}


Job::JobResult EDFT::EvalJob(Topology *top, Job *job, QMThread *opThread) {

    string output;
    
    bool _run_status;
    bool _parse_log_status;
    bool _parse_orbitals_status;

   
    
    Orbitals _orbitals;
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
    XTP_LOG(logINFO,*pLog) << TimeStamp() << " Evaluating site " << seg->getId() << flush; 

    // log, com, and orbital files will be stored in ORB_FILES/package_name/frame_x/mol_ID/
    // extracted information will be stored in  ORB_FILES/molecules/frame_x/molecule_ID.orb
    
    string edft_work_dir = "OR_FILES";
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());      
    string ID   = boost::lexical_cast<string>( seg->getId() );

    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
    
   _qmpackage->setLog( pLog );  
   _qmpackage->Initialize( &_package_options );



   
      string   orbitals_storage_dir = "molecules";
      string   qmpackage_work_dir  = edft_work_dir + "/" + _package + "/" + frame_dir + "/mol_" + ID;
 

    string ORB_DIR = edft_work_dir + "/" + orbitals_storage_dir + "/" + frame_dir;
    boost::filesystem::create_directories( qmpackage_work_dir );     
    boost::filesystem::create_directories(ORB_DIR); 
    
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
            XTP_LOG(logERROR,*pLog) << _package << " run failed" << flush;
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
            XTP_LOG(logERROR,*pLog) << "QM log incomplete" << flush;
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
            XTP_LOG(logERROR,*pLog) << "QM orbitals not parsed" << flush;
            jres.setOutput( output ); 
            jres.setStatus(Job::FAILED);
            delete _qmpackage;
            return jres;
        } else {
            output += "orbitals parsed; " ;
        }
    }
 
   // Trim virtual orbitals
    if ( _do_trim ) {
        
       int factor = 2;
       if ( !_do_parse ) { // orbitals must be loaded from a file
           boost::filesystem::path arg_path;
           string ORB_FILE = ( arg_path / ORB_DIR / (format("molecule_%1%.orb") % ID ).str() ).c_str() ;
           XTP_LOG(logDEBUG,*pLog) << "Loading orbitals from " << ORB_FILE << flush;  
           try{
               _orbitals.ReadFromCpt(ORB_FILE);
           }
           catch(std::runtime_error& error){
               XTP_LOG(logERROR,*pLog) << "Failed loading orbitals from " << ORB_FILE << flush; 
               output += "failed loading " + ORB_FILE;
               jres.setOutput( output ); 
               jres.setStatus(Job::FAILED);
               delete _qmpackage;
               return jres;
           }
           
        }        
       
       _orbitals.Trim(factor);   
        XTP_LOG(logDEBUG,*pLog) << "Trimming virtual orbitals from " 
         << _orbitals.getNumberOfLevels() - _orbitals.getNumberOfElectrons() << " to " 
         << _orbitals.getNumberOfElectrons()*factor << flush;   
       output += "orbitals trimmed; " ;
    }   

   
   if ( _do_parse ){
    // save orbitals
    
    string ORB_FILE = "molecule_" + ID + ".orb";
    XTP_LOG(logDEBUG,*pLog) << "Serializing to " <<  ORB_FILE << flush;
    _orbitals.WriteToCpt(ORB_DIR+"/"+ORB_FILE);
    // ofs.close();
    
     if(_qmpackage->getPackageName()=="orca"){
            XTP_LOG(logINFO,*pLog) << "Copying monomer .gbw file to orb folder" << flush;
            string   qmpackage_gbw_dir  = edft_work_dir + "/" + _package + "/" + frame_dir + "/mol_" + ID+"/system.gbw";           
            string gbwFile  = ORB_DIR+"/"+(format("%1%_%2%%3%") % "molecule" % ID % ".gbw").str();
            boost::filesystem::copy_file(qmpackage_gbw_dir, gbwFile,boost::filesystem::copy_option::overwrite_if_exists);
    }
    
    XTP_LOG(logDEBUG,*pLog) << "Done serializing " <<  ORB_FILE << flush;
   }
   
  
   
   // Clean run
   _qmpackage->CleanUp();
   delete _qmpackage;
        
    XTP_LOG(logINFO,*pLog) << TimeStamp() << " Finished evaluating site " << seg->getId() << flush; 
 
    Property _job_summary;
        Property *_output_summary = &_job_summary.add("output","");
        Property *_segment_summary = &_output_summary->add("segment","");
         string segName = seg->getName();
         segId = seg->getId();
        _segment_summary->setAttribute("id", segId);
        _segment_summary->setAttribute("type", segName);
        _segment_summary->setAttribute("homo", _orbitals.getEnergy( _orbitals.getNumberOfElectrons() ));
        _segment_summary->setAttribute("lumo", _orbitals.getEnergy( _orbitals.getNumberOfElectrons() + 1 ));
    
    // output of the JOB 
    jres.setOutput( _job_summary );
    jres.setStatus(Job::COMPLETE);

    // dump the LOG
    //cout << *pLog;
    
    return jres;

}

}}

#endif	/* _VOTCA_XTP_EDFT_H */
