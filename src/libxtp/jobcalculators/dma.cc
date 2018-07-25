/*
 *            Copyright 2009-2017 The VOTCA Development Team
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


#include "dma.h"

#include <boost/algorithm/string/replace.hpp>
#include <votca/xtp/qminterface.h>


using boost::format;

namespace votca { namespace xtp {

/**
* \brief Distributed multipole analysis using Gaussian input
*
* Evaluates distributed multipoles
* Requires GAUSSIAN and GDMA
*
* Callname: dma
*/

void DMA::Initialize(Property *options) {

    _do_input = false;
    _do_orbitals = false;
    _do_dma = false;
           
    _maverick = (_nThreads == 1) ? true : false;
    
    string key = "options." + Identify();
    string _package_xml = options->get(key+".package").as<string> ();
    _chkFile = options->get(key+".chk").as<string> ();
    _executable = options->get(key+".executable").as<string> ();

    _density = options->get(key+".density").as<string> ();
    _limit = options->get(key+".multipoles.limit").as<int> ();
    _radius = options->get(key+".multipoles.radius").as<double> ();
    _switch = options->get(key+".multipoles.switch").as<int> ();
    _outFile = options->get(key+".output").as<string> ();
    
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

void DMA::WriteJobFile(ctp::Topology *top) {

    cout << endl << "... ... Writing job file: " << flush;
    ofstream ofs;
    ofs.open(_jobfile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);
 
    ofs << "<jobs>" << endl;   

    std::vector<ctp::Segment*> &segments = top->Segments();
    std::vector<ctp::Segment*>::iterator sit;
 
    int jobCount = 0;

    for (sit = segments.begin(); sit != segments.end(); ++sit) {
    
        int id = ++jobCount;
        string tag = "";

        Property Input;
        Property *pInput = &Input.add("input","");
        
        int segment_id =  (*sit)->getId();
        Property *pSegment =  &pInput->add("segment" , boost::lexical_cast<string>(segment_id)  ) ;
        pSegment->setAttribute<string>("type", (*sit)->getName() );
        pSegment->setAttribute<int>("id", segment_id );
        ctp::Job job(id, tag, Input, ctp::Job::AVAILABLE );
        job.ToStream(ofs,"xml");
    }
     

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    
    cout << jobCount << " jobs" << flush;
    
}


ctp::Job::JobResult DMA::EvalJob(ctp::Topology *top, ctp::Job *job, ctp::QMThread *opThread) {

    string output;
    
    bool _orbitals_status;
    //bool _formchk_status;
    bool _dma_status;

    ctp::Job::JobResult jres = ctp::Job::JobResult();
    Property _job_input = job->getInput();  
    list<Property*> lSegments = _job_input.Select( "segment" );  
    
    int segId = lSegments.front()->getAttribute<int>( "id" );
    string segType = lSegments.front()->getAttribute<string>( "type" );
    
    ctp::Segment *seg = top->getSegment( segId );
    assert( seg->getName() == segType ); 
    std::vector<ctp::Segment*> segments;
    segments.push_back( seg );
    QMInterface interface;
    Orbitals orbital;
    orbital.QMAtoms()=interface.Convert(segments);
    ctp::Logger* pLog = opThread->getLogger();
    CTP_LOG(ctp::logINFO,*pLog) << ctp::TimeStamp() << " Evaluating site " << seg->getId() << flush; 

    // log, com, and orbital files will be stored in ORB_FILES/package_name/frame_x/mol_ID/
    // extracted information will be stored in  ORB_FILES/molecules/frame_x/molecule_ID.orb
    
    string dma_work_dir = "DMA_FILES";
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());      
    string ID   = boost::lexical_cast<string>( seg->getId() );

    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );
    
   _qmpackage->setLog( pLog );  
   _qmpackage->Initialize( _package_options );

   string qmpackage_work_dir;
   
   qmpackage_work_dir  = dma_work_dir + "/" + _package + "/" + frame_dir + "/mol_" + ID;
   boost::filesystem::create_directories( qmpackage_work_dir );     
    
   _qmpackage->setRunDir( qmpackage_work_dir );
   
    // if asked, prepare the input files
    if ( _do_input ) {
        _qmpackage->WriteInputFile( orbital );
    }
        
   // Run the Gaussian executable
    if ( _do_orbitals ) {
        _orbitals_status = _qmpackage->Run( orbital);
        if ( !_orbitals_status ) {
            output += "run failed; " ;
            CTP_LOG(ctp::logERROR,*pLog) << _package << " run failed" << flush;
            jres.setOutput( output ); 
            jres.setStatus(ctp::Job::FAILED);
            delete _qmpackage;
            return jres;
        } else {
            output += "orbitals completed; " ;
        }
    }
        
   // Parse log files
    if ( _do_dma ) {
        
        
        CTP_LOG(ctp::logDEBUG,*pLog) << "DMA: running formcheck on [" << _chkFile << "]" << flush;

        if (std::system(NULL)) {
            string _command;
            _command  = "cd " + qmpackage_work_dir + "; formchk " + _chkFile;

            int i = std::system ( _command.c_str() );
 
            // prepare the GDMA input file
            //Title ""
            // File AAE007_neutral.fchk
                //Multipoles
                // Limit 2
                // Limit 1 H
                //  Radius H 0.35
                //  Punch workout_neutral_6-311+gdp.mps
            //Start
            //Finish
            
            std::fstream fs;
            string _dma_input_file = qmpackage_work_dir + "/gdma.in";
            
            cout << _dma_input_file << std::endl;
            
            fs.open (_dma_input_file.c_str(), std::fstream::out );
            if (fs.is_open()) {
                fs << "Title " << "GDMA\n";
                fs << "Density " << _density << std::endl;

                // chk -> fchk
                string _fchkFile = _chkFile;
                boost::algorithm::replace_all(_fchkFile, ".chk", ".fchk");  

                fs << "File " << _fchkFile << std::endl;
                fs << "Multipoles\n"  ;
                fs << "Switch " << _switch << std::endl;
                fs << "Limit "  << _limit << std::endl;
                fs << "Radius " << _radius << std::endl;
                fs << "Punch " <<  _outFile << std::endl ;            
                fs << "Start\n"  ;
                fs << "Finish\n"  ;
                fs.close();            

                 _command  = "cd " + qmpackage_work_dir + "; " + _executable + " < gdma.in > gdma.out";

                CTP_LOG(ctp::logINFO,*pLog) << _command << flush;

                int j = std::system ( _command.c_str() );

                _dma_status = (i == 0) && (j == 0);
            
            }
            else {
                CTP_LOG(ctp::logDEBUG,*pLog) << "Error opening file " << _dma_input_file << flush;
                _dma_status = false;
            }

            if ( !_dma_status ) {
                output += "DMA job incomplete; ";
                CTP_LOG(ctp::logERROR,*pLog) << "DMA job incomplete" << flush;
                jres.setOutput( output ); 
                jres.setStatus(ctp::Job::FAILED);
               return jres;
            } else {
                output += "DMA complete; " ;
            }
        }
        
 
    }

   
   // Clean run
   _qmpackage->CleanUp();
   delete _qmpackage;
        
    CTP_LOG(ctp::logINFO,*pLog) << ctp::TimeStamp() << " Finished evaluating site " << seg->getId() << flush; 
 
    Property _job_summary;
        Property *_output_summary = &_job_summary.add("output","");
        Property *_segment_summary = &_output_summary->add("segment","");
         string segName = seg->getName();
         segId = seg->getId();
        _segment_summary->setAttribute("id", segId);
    
    // output of the JOB 
    jres.setOutput( _job_summary );
    jres.setStatus(ctp::Job::COMPLETE);

    
    return jres;

}

}}
