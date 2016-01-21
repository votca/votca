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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/xtp/votca_xtp_config.h>

#include "egwbse.h"



#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
// #include <votca/xtp/logger.h>
#include <votca/xtp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        // +++++++++++++++++++++++++++++ //
        // GWBSE MEMBER FUNCTIONS         //
        // +++++++++++++++++++++++++++++ //

        void EGWBSE::CleanUp() {

        }

        void EGWBSE::Initialize(Property *options) {


            // tasks to be done by IBSE: dft_input, dft_run, dft_parse, mgbft, bse_coupling
            _do_dft_input = false;
            _do_dft_run = false;
            _do_dft_parse = false;
            _do_gwbse = false;

            // update options with the VOTCASHARE defaults   
            UpdateWithDefaults(options);
            ParseOptionsXML(options);

            // register all QM packages (Gaussian, turbomole, etc))
            QMPackageFactory::RegisterAll();



        }

        void EGWBSE::ParseOptionsXML(Property* options) {



            _maverick = (_nThreads == 1) ? true : false;
            string key = "options." + Identify();
            // job tasks
            string _tasks_string = options->get(key + ".tasks").as<string> ();
            if (_tasks_string.find("input") != std::string::npos) _do_dft_input = true;
            if (_tasks_string.find("dft") != std::string::npos) _do_dft_run = true;
            if (_tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
            if (_tasks_string.find("gwbse") != std::string::npos) _do_gwbse = true;


            key = "options." + Identify() + ".job";
            _jobfile = options->get(key + ".file").as<string>();

            // options for gwbse
            key = "options." + Identify();
            string _gwbse_xml = options->get(key + ".gwbse").as<string> ();
            load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());



            // options for dft package
            string _package_xml = options->get(key + ".package").as<string> ();
            //cout << endl << "... ... Parsing " << _package_xml << endl ;
            load_property_from_xml(_package_options, _package_xml.c_str());
            key = "package";
            _package = _package_options.get(key + ".name").as<string> ();


        }
        
        
        void EGWBSE::WriteJobFile(Topology *top) {

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

        Job::JobResult EGWBSE::EvalJob(Topology *top, Job *job, QMThread *opThread) {

            cout << "Starting GW-BSE";
            Orbitals _orbitals;
            Job::JobResult jres = Job::JobResult();
            Property _job_input = job->getInput();
            list<Property*> lSegments = _job_input.Select("segment");

            vector < Segment* > segments;
            int segId = lSegments.front()->getAttribute<int>("id");
            string segType = lSegments.front()->getAttribute<string>("type");

            Segment *seg = top->getSegment(segId);
            assert(seg->Name() == segType);
            segments.push_back(seg);

            Logger* pLog = opThread->getLogger();
            LOG(logINFO, *pLog) << TimeStamp() << " Evaluating site " << seg->getId() << flush;

            
            string output;
            
            // directories and files
            path arg_path;
            string egwbse_work_dir = "OR_FILES";
            string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());     
            string orb_file = (format("%1%_%2%%3%") % "molecule" % segId % ".orb").str();
            string _mol_dir = ( format("%1%%2%%3%") % "molecule" % "_" % segId ).str();
            string _package_append = _package + "_gwbse";
            string _qmpackage_work_dir = (arg_path / egwbse_work_dir / _package_append / frame_dir / _mol_dir).c_str();

    
            // get the corresponding object from the QMPackageFactory
            QMPackage *_qmpackage = QMPackages().Create(_package);
            // set a log file for the package
            _qmpackage->setLog(pLog);
            // set the run dir 
            _qmpackage->setRunDir(_qmpackage_work_dir);
            // get the package options
            _qmpackage->Initialize(&_package_options);
            
            
            // different tasks
            
            // create input for DFT
            if ( _do_dft_input ){
                boost::filesystem::create_directories( _qmpackage_work_dir );
                _qmpackage->WriteInputFile( segments );
            }
            
            bool _run_dft_status = false;
            if (_do_dft_run) {
                _run_dft_status = _qmpackage->Run();
                if (!_run_dft_status) {
                    output += "run failed; ";
                    LOG(logERROR, *pLog) << _qmpackage->getPackageName() << " run failed" << flush;
                    cout << *pLog;
                    jres.setOutput(output);
                    jres.setStatus(Job::FAILED);
                    delete _qmpackage;
                    return jres;
                }
            }

            // parse the log/orbitals files
            bool _parse_log_status = false;
            bool _parse_orbitals_status = false;
            if (_do_dft_parse) {
                _parse_log_status = _qmpackage->ParseLogFile(&_orbitals);

                if (!_parse_log_status) {
                    output += "log incomplete; ";
                    LOG(logERROR, *pLog) << "LOG parsing failed" << flush;
                    cout << *pLog;
                    jres.setOutput(output);
                    jres.setStatus(Job::FAILED);
                    delete _qmpackage;
                    return jres;
                }

                _parse_orbitals_status = _qmpackage->ParseOrbitalsFile(&_orbitals);

                if (!_parse_orbitals_status) {
                    output += "orbfile failed; ";
                    LOG(logERROR, *pLog) << "Orbitals parsing failed" << flush;
                    cout << *pLog;
                    jres.setOutput(output);
                    jres.setStatus(Job::FAILED);
                    delete _qmpackage;
                    return jres;
                }
            } // end of the parse orbitals/log


            
 
            if (_do_gwbse) {

                if (!_do_dft_parse) {

                    // load the DFT data from serialized orbitals object
                    string DIR = egwbse_work_dir + "/molecules_gwbse/" + frame_dir;
                    std::ifstream ifs((DIR + "/" + orb_file).c_str());
                    LOG(logDEBUG, *pLog) << TimeStamp() << " Loading DFT data from " << DIR << "/" << orb_file << flush;
                    boost::archive::binary_iarchive ia(ifs);
                    ia >> _orbitals;
                    ifs.close();
                }
                GWBSE _gwbse; 
                _gwbse.Initialize(&_gwbse_options);
                // _gwbse.setLogger(pLog);
                
                
                // define own logger for GW-BSE that is written into a runFolder logfile
                Logger gwbse_logger(logDEBUG);
                gwbse_logger.setMultithreading(false);
                _gwbse.setLogger(&gwbse_logger);
                gwbse_logger.setPreface(logINFO,    (format("\nGWBSE INF ...") ).str());
                gwbse_logger.setPreface(logERROR,   (format("\nGWBSE ERR ...") ).str());
                gwbse_logger.setPreface(logWARNING, (format("\nGWBSE WAR ...") ).str());
                gwbse_logger.setPreface(logDEBUG,   (format("\nGWBSE DBG ...") ).str());
                
                
                //bool _evaluate = _gwbse.Evaluate(&_orbitals);
                _gwbse.Evaluate(&_orbitals);
                
                // write logger to log file
                ofstream ofs;
                string gwbse_logfile = _qmpackage_work_dir + "/gwbse.log";
                ofs.open(gwbse_logfile.c_str(), ofstream::out);
                if (!ofs.is_open()) {
                    throw runtime_error("Bad file handle: " + gwbse_logfile);
                }    
                ofs << gwbse_logger << endl;
                ofs.close();

            }

            LOG(logINFO, *pLog) << TimeStamp() << " Finished evaluating site " << seg->getId() << flush;


            
            
            LOG(logDEBUG, *pLog) << "Saving data to " << orb_file << flush;
            string DIR = egwbse_work_dir + "/molecules_gwbse/" + frame_dir;
            boost::filesystem::create_directories(DIR);  
            std::ofstream ofs((DIR + "/" + orb_file).c_str());
            boost::archive::binary_oarchive oa(ofs);

            //  if ( !( _store_orbitals && _do_parse && _parse_orbitals_status) )   _store_orbitals = false;
            //  if ( !( _store_overlap && _do_parse && _parse_log_status) )    _store_overlap = false;
            //  if ( !( _store_integrals && _do_project && _calculate_integrals) )  {
            //      _store_integrals = false; 
            //  } else {
            //      _orbitalsAB.setIntegrals( &_JAB );
            //  }

            //  _orbitalsAB.setStorage( _store_orbitals, _store_overlap, _store_integrals );


            oa << _orbitals;
            ofs.close();





            Property _job_summary;
            Property *_output_summary = &_job_summary.add("output", "");
            Property *_segment_summary = &_output_summary->add("segment", "");
            string segName = seg->getName();
            segId = seg->getId();
            _segment_summary->setAttribute("id", segId);
            _segment_summary->setAttribute("type", segName);
            // output of the JOB 
            jres.setOutput(_job_summary);
            jres.setStatus(Job::COMPLETE);

            // dump the LOG

            return jres;
        }

    }


};
