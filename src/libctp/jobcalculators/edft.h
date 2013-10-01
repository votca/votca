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


#ifndef _VOTCA_CTP_EDFT_H
#define	_VOTCA_CTP_EDFT_H

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
#include <boost/format.hpp>

using boost::format;

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

    string   Identify() { return "edft"; }
    void    Initialize(Property *options);
    void    ParseOrbitalsXML(Topology *top, Property *options);
    void    WriteJobList(Topology *top);
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *thread);

    void    CleanUp();


private:

    // what to do
    bool                _do_input;
    bool                _do_run;
    bool                _do_parse;
    bool                _do_trim;
    // conversion to GW (might go away from edft again)
    bool                _do_convert;
    bool                _do_gwbse_input;
    bool                _do_gwbse_run;
    bool                _do_gwbse_parse;
    // what to write in the storage
    bool                _store_orbitals;
    bool                _store_qppert;
    bool                _store_qpdiag;
    bool                _store_singlets;
    bool                _store_triplets;

    //bool   _maverick;
    string _outParent;
    string _outMonDir;
    
    string _package;
    Property _package_options;   
    
    string _gwpackage;
    Property _gwpackage_options;   

};

void EDFT::CleanUp() {

}

void EDFT::Initialize(Property *options) {

    _do_input = false;
    _do_run = false;
    _do_parse = false;
    _do_trim = false;
    
    // conversion to GW 
    _do_convert = false;
    _do_gwbse_input = false;
    _do_gwbse_run = false;
    _do_gwbse_parse = false;
        
    _maverick = (_nThreads == 1) ? true : false;
    
    string key = "options." + Identify();
    _jobfile = options->get(key + ".job_file").as<string>();

    string _package_xml = options->get(key+".package").as<string> ();

    
    string _tasks_string = options->get(key+".tasks").as<string> ();
    if (_tasks_string.find("input") != std::string::npos) _do_input = true;
    if (_tasks_string.find("run") != std::string::npos) _do_run = true;
    if (_tasks_string.find("trim") != std::string::npos) _do_trim = true;
    if (_tasks_string.find("parse") != std::string::npos) _do_parse = true;    
    // GW-BSE tasks
    if (_tasks_string.find("convert") != std::string::npos) _do_convert = true;   
    if (_tasks_string.find("gwbse_setup") != std::string::npos) _do_gwbse_input = true;
    if (_tasks_string.find("gwbse_exec") != std::string::npos) _do_gwbse_run = true;    
    if (_tasks_string.find("gwbse_read") != std::string::npos) _do_gwbse_parse = true;
    
    string _store_string = options->get(key+".store").as<string> ();
    if (_store_string.find("orbitals") != std::string::npos) _store_orbitals = true;
    if (_store_string.find("qppert") != std::string::npos) _store_qppert = true;
    if (_store_string.find("qpdiag") != std::string::npos) _store_qpdiag = true;
    if (_store_string.find("singlets") != std::string::npos) _store_singlets = true;
    if (_store_string.find("triplets") != std::string::npos) _store_triplets = true;
    
    load_property_from_xml( _package_options, _package_xml.c_str() );    
    key = "package";
    _package = _package_options.get(key+".name").as<string> ();


   
    // only required, if GWBSE is to be run
    if ( _do_gwbse_input || _do_gwbse_run || _do_gwbse_parse ){
        key = "options." + Identify();
        string _gwpackage_xml = options->get(key+".gwpackage").as<string> ();
        load_property_from_xml( _gwpackage_options, _gwpackage_xml.c_str() );  
        key = "package";
        _gwpackage = _gwpackage_options.get(key+".name").as<string> ();
    }
    
    
    // register all QM packages (Gaussian, turbomole, nwchem))
    QMPackageFactory::RegisterAll(); 

}

void EDFT::WriteJobList(Topology *top) {

    string jobFile = "edft.jobs";   
    
    ofstream ofs;
    ofs.open(jobFile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("Bad file handle: " + jobFile);
 
    ofs << "<jobs>" << endl;   

    /* this is only good when ALL molecules shall be written out 
    for (vector<Segment*>::iterator sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {
        int id = (*sit)->getId();
        string tag = (*sit)->getId();
        string input = "";
        string stat = "AVAILABLE";
        Job job(id, tag, input, stat);
        job.ToStream(ofs,"xml");
    }
    */

    QMNBList::iterator pit;
    QMNBList &nblist = top->NBList();    

    int jobCount = 0;
    if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
    } 

    map< int,Segment* > segments;
    map< int,Segment* >::iterator sit;

    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {
        
        int id1 = (*pit)->Seg1()->getId();
        int id2 = (*pit)->Seg2()->getId();
	segments[id1] = (*pit)->Seg1();
        segments[id2] = (*pit)->Seg2();

    }

    for (sit = segments.begin(); sit != segments.end(); ++sit) {
    
        int id = ++jobCount;
        
        string tag = (format("%1$s") % sit->first).str();
        string input = sit->second->getName();
        string stat = "AVAILABLE";
        Job job(id, tag, input, stat);
        job.ToStream(ofs,"xml");
    }
     

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    
}


Job::JobResult EDFT::EvalJob(Topology *top, Job *job, QMThread *opThread) {

    string output;
    
    bool _run_status;
    bool _parse_log_status;
    bool _parse_orbitals_status;
    bool _convert_status;


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
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());      
    string ID   = boost::lexical_cast<string>( seg->getId() );

    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
    
   _qmpackage->setLog( pLog );  
   _qmpackage->Initialize( &_package_options );

   string orbitals_storage_dir;
   string qmpackage_work_dir;
   
   /* To avoid possible mixup between standard DFT runs, e.g., for coupling
    * elements and GWBSE DFT basis runs, modify storage files/directories */
   if ( _qmpackage->ECPRequested() ){
       orbitals_storage_dir = "molecules_gwbse";
       qmpackage_work_dir  = edft_work_dir + "/" + _package + "_gwbse/" + frame_dir + "/mol_" + ID;
   } else {
       orbitals_storage_dir = "molecules";
       qmpackage_work_dir  = edft_work_dir + "/" + _package + "/" + frame_dir + "/mol_" + ID;
   }

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
 
   // Trim virtual orbitals
    if ( _do_trim ) {
        
       int factor = 2;
       if ( !_do_parse ) { // orbitals must be loaded from a file
           string ORB_FILE = "molecule_" + ID + ".orb";
           LOG(logDEBUG,*pLog) << "Loading orbitals from " << ORB_FILE << flush;    
           std::ifstream ifs( (ORB_DIR + "/" + ORB_FILE).c_str() );
           boost::archive::binary_iarchive ia( ifs );
           ia >> _orbitals;
           ifs.close();
       }        
       
       _orbitals.Trim(factor);   
        LOG(logDEBUG,*pLog) << "Trimming virtual orbitals from " 
         << _orbitals.getNumberOfLevels() - _orbitals.getNumberOfElectrons() << " to " 
         << _orbitals.getNumberOfElectrons()*factor << flush;   
       output += "orbitals trimmed; " ;
    }   

   
   if ( _do_parse ){
    // save orbitals
    string ORB_FILE = "molecule_" + ID + ".orb";
    LOG(logDEBUG,*pLog) << "Serializing to " <<  ORB_FILE << flush;
    std::ofstream ofs( (ORB_DIR + "/" + ORB_FILE).c_str() );
    boost::archive::binary_oarchive oa( ofs );
    oa << _orbitals;
    ofs.close();
   }
   
   // Convert to GW
    if ( _do_convert ) {

        if ( !_do_parse ) { // orbitals must be loaded from a file
           string ORB_FILE = "molecule_" + ID + ".orb";
           LOG(logDEBUG,*pLog) << "Loading orbitals from " << ORB_FILE << flush;    
           std::ifstream ifs( (ORB_DIR + "/" + ORB_FILE).c_str() );
           boost::archive::binary_iarchive ia( ifs );
           ia >> _orbitals;
           ifs.close();
       } 
        
       _convert_status = _qmpackage->ConvertToGW( &_orbitals );
       if ( !_convert_status ) {
            output += "conversion failed; " ;
            LOG(logERROR,*pLog) << _package << " conversion failed" << flush;
            jres.setOutput( output ); 
            jres.setStatus(Job::FAILED);
            delete _qmpackage;
            return jres;
        } else {
            output += "conversion completed; " ;
        }
        
    }
   

            // to run any of the GW-BSE steps, initialize QMPackage
            if (_do_gwbse_input || _do_gwbse_run || _do_gwbse_parse) {

                // initialize GW as another QMPackage
                QMPackage *_gwpackage = QMPackages().Create("gw");
                _gwpackage->Initialize(&_gwpackage_options);
                _gwpackage->setLog(pLog);
                _gwpackage->setRunDir(qmpackage_work_dir);


                if (_do_gwbse_input || _do_gwbse_parse ) {

                    // DFT data must be read from storage or have been parsed
                    if (!_do_parse) { // orbitals must be loaded from a file
                        string ORB_FILE = "molecule_" + ID + ".orb";
                        LOG(logDEBUG, *pLog) << "Loading orbitals from " << ORB_FILE << flush;
                        std::ifstream ifs((ORB_DIR + "/" + ORB_FILE).c_str());
                        boost::archive::binary_iarchive ia(ifs);
                        ia >> _orbitals;
                        ifs.close();
                    }

                    _gwpackage->WriteInputFile(segments, &_orbitals);

                }

                if (_do_gwbse_run) {
                    // run the GWBSE executable
                    _run_status = _gwpackage->Run();
                    if (!_run_status) {
                        output += "run failed; ";
                        LOG(logERROR, *pLog) << " GWBSE run failed" << flush;
                        jres.setOutput(output);
                        jres.setStatus(Job::FAILED);
                        delete _gwpackage;
                        return jres;
                    } else {
                        output += "run completed; ";
                    }
                }

                if (_do_gwbse_parse) {
                    // parse and store output from GWBSE
                    _parse_log_status = _gwpackage->ParseLogFile(&_orbitals);
                    if (!_parse_log_status) {
                        output += "log incomplete; ";
                        LOG(logERROR, *pLog) << "GWBSE log incomplete" << flush;
                        jres.setOutput(output);
                        jres.setStatus(Job::FAILED);
                        delete _gwpackage;
                        return jres;
                    } else {
                        output += "log parsed; ";
                        // save orbitals
                        string ORB_FILE = "molecule_" + ID + ".orb";
                        LOG(logDEBUG,*pLog) << "Serializing to " <<  ORB_FILE << flush;
                        std::ofstream ofs( (ORB_DIR + "/" + ORB_FILE).c_str() );
                        boost::archive::binary_oarchive oa( ofs );
                        oa << _orbitals;
                        ofs.close();
                    }

             
                }
            }
   
   // Clean run
   _qmpackage->CleanUp();
   delete _qmpackage;
        
    // GENERATE OUTPUT AND FORWARD TO PROGRESS OBSERVER (RETURN)
    jres.setStatus(Job::COMPLETE);


                                
    LOG(logINFO,*pLog) << TimeStamp() << " Finished evaluating site " << seg->getId() << flush; 
    
    // output of the JOB 
    jres.setOutput( output );

    // dump the LOG
    cout << *pLog;
    
    return jres;
    
   
}

}}

#endif	/* _VOTCA_CTP_EDFT_H */
