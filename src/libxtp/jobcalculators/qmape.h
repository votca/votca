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

#ifndef VOTCA_XTP_QMAPECALC_H
#define	VOTCA_XTP_QMAPECALC_H

#include <votca/xtp/pewald3d.h>
#include <votca/xtp/parallelxjobcalc.h>
#include <votca/xtp/xmapper.h>
#include <votca/xtp/xjob.h>
#include <votca/xtp/xinductor.h>
#include <votca/xtp/xinteractor.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmapemachine.h>
#include <boost/format.hpp>

using boost::format;

namespace votca { namespace xtp {

   
class QMAPE : public xtp::ParallelXJobCalc< vector<xtp::Job*>, xtp::Job*, xtp::Job::JobResult >
{

public:

    QMAPE() {};
   ~QMAPE() {};
   
    string          Identify() { return "qmape"; }
    void            Initialize(Property *);

    void            CustomizeLogger(xtp::QMThread *thread);
    void            PreProcess(xtp::Topology *top);
    xtp::Job::JobResult  EvalJob(xtp::Topology *top, xtp::Job *job, xtp::QMThread *thread);
    xtp::XJob            ProcessInputString(xtp::Job *job, xtp::Topology *top, xtp::QMThread *thread);

private:
    
    // ======================================== //
    // MULTIPOLE ALLOCATION, XJOBS, ADD. OUTPUT //
    // ======================================== //

	string                         _xml_file;
	string                         _mps_table;
	string                         _polar_bg_arch;
	xtp::XMpsMap                   _mps_mapper;
	bool                           _pdb_check;
	bool                           _ptop_check;
    

    // Multipole Interaction parameters
    string                          _method;
    Property                        _dft_opt;
    // GWBSE options
    string                          _gwbse;
    Property                        _gwbse_opt;
    // XJob logbook (file output)
    string                          _outFile;
   
    
    Property                       *_options;
};

// ========================================================================== //
//                      PARALLELCALC MEMBER FUNCTIONS                         //
// ========================================================================== //


void QMAPE::Initialize(Property *options) {
    
	_options = options;

    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;

    _maverick = (_nThreads == 1) ? true : false;

    string key = "options."+Identify();

        
    _jobfile = options->ifExistsReturnElseThrowRuntimeError<string>(key+".job_file");
        

	key = "options."+Identify()+".ewald";	
        _xml_file = options->ifExistsReturnElseThrowRuntimeError<string>(key+".mapping");
        _mps_table = options->ifExistsReturnElseThrowRuntimeError<string>(key+".mps_table");
        _polar_bg_arch = options->ifExistsReturnElseReturnDefault<string>(key+".polar_bg","");
        _pdb_check = options->ifExistsReturnElseReturnDefault<bool>(key+".pdb_check",false);
	_ptop_check = options->ifExistsReturnElseReturnDefault<bool>(key+".ptop_check",false);

    
    key = "options."+Identify();
    
      
    string package_xml = options->ifExistsReturnElseThrowRuntimeError<string>(key+".dft.dftengine");
    load_property_from_xml(_dft_opt, package_xml.c_str());
       
    key = "options."+Identify()+".gwbse";
    if ( options->exists(key)) { 
    	cout << endl << "... ... Configure for excited states (DFT+GWBSE)" << flush;
        string gwbse_xml = options->ifExistsReturnElseThrowRuntimeError<string>(key+".gwbse_options");
        load_property_from_xml(_gwbse_opt, gwbse_xml.c_str());
    }
    else {
        cout << endl << "... ... Configure for ground states (DFT)" << flush;
    }
    return;
}


void QMAPE::PreProcess(xtp::Topology *top) {
    // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
    cout << endl << "... ... Initialize MPS-mapper: " << flush;
    _mps_mapper.GenerateMap(_xml_file, _mps_table, top);
    return;
}


void QMAPE::CustomizeLogger(xtp::QMThread *thread) {
    
    // CONFIGURE LOGGER
    xtp::Logger* log = thread->getLogger();
    log->setReportLevel(xtp::logDEBUG);
    log->setMultithreading(_maverick);

    log->setPreface(xtp::logINFO,    (format("\nT%1$02d INF ...") % thread->getId()).str());
    log->setPreface(xtp::logERROR,   (format("\nT%1$02d ERR ...") % thread->getId()).str());
    log->setPreface(xtp::logWARNING, (format("\nT%1$02d WAR ...") % thread->getId()).str());
    log->setPreface(xtp::logDEBUG,   (format("\nT%1$02d DBG ...") % thread->getId()).str()); 
    return;
}


// ========================================================================== //
//                            QMAPE MEMBER FUNCTIONS                          //
// ========================================================================== //


xtp::XJob QMAPE::ProcessInputString(xtp::Job *job,xtp::Topology *top, xtp::QMThread *thread) {

    // Input string looks like this:
    // <id1>:<name1>:<mpsfile1> <id2>:<name2>: ... ... ...

    string input = job->getInput().as<string>();
    vector<xtp::Segment*> qmSegs;
    vector<string>   qmSegMps;
    vector<string> split;
    Tokenizer toker(input, " \t\n");
    toker.ToVector(split);

    for (unsigned i = 0; i < split.size(); ++i) {

        string id_seg_mps = split[i];
        vector<string> split_id_seg_mps;
        Tokenizer toker(id_seg_mps, ":");
        toker.ToVector(split_id_seg_mps);

        int segId = boost::lexical_cast<int>(split_id_seg_mps[0]);
        string segName = split_id_seg_mps[1];
        string mpsFile = split_id_seg_mps[2];

        xtp::Segment *seg = top->getSegment(segId);
        if (seg->getName() != segName) {
            XTP_LOG(xtp::logERROR,*(thread->getLogger()))
                << "ERROR: Seg " << segId << ":" << seg->getName() << " "
                << " maltagged as " << segName << ". Skip job ..." << flush;
            throw std::runtime_error("Input does not match topology.");
        }

        qmSegs.push_back(seg);
        qmSegMps.push_back(mpsFile);
    }

    return xtp::XJob(job->getId(), job->getTag(), qmSegs, qmSegMps, top);
}


xtp::Job::JobResult QMAPE::EvalJob(xtp::Topology *top, xtp::Job *job, xtp::QMThread *thread) {
    
    // SILENT LOGGER FOR QMPACKAGE
    xtp::Logger* log = thread->getLogger();    
    xtp::Logger qlog;
    qlog.setReportLevel(xtp::logDEBUG);
    qlog.setMultithreading(_maverick);
    qlog.setPreface(xtp::logINFO,    (format("\nQ%1$02d ... ...") % thread->getId()).str());
    qlog.setPreface(xtp::logERROR,   (format("\nQ%1$02d ERR ...") % thread->getId()).str());
    qlog.setPreface(xtp::logWARNING, (format("\nQ%1$02d WAR ...") % thread->getId()).str());
    qlog.setPreface(xtp::logDEBUG,   (format("\nQ%1$02d DBG ...") % thread->getId()).str());

    // CREATE XJOB FROM JOB INPUT STRING
    XTP_LOG(xtp::logINFO,*log)
        << "Job input = " << job->getInput().as<string>() << flush;
    xtp::XJob xjob = this->ProcessInputString(job, top, thread);  

	// SETUP POLAR TOPOLOGY (GENERATE VS LOAD IF PREPOLARIZED)
	if (_polar_bg_arch == "") {
		XTP_LOG(xtp::logINFO,*log) << "Mps-Mapper: Generate FGC FGN BGN" << flush;
		_mps_mapper.Gen_FGC_FGN_BGN(top, &xjob, thread);
	}
	else {
		XTP_LOG(xtp::logINFO,*log) << "Mps-Mapper: Generate FGC, load FGN BGN from '"
				<< _polar_bg_arch << "'" << flush;
		_mps_mapper.Gen_FGC_Load_FGN_BGN(top, &xjob, _polar_bg_arch, thread);
	}
    XTP_LOG(xtp::logINFO,*log) << xjob.getPolarTop()->ShellInfoStr() << flush;

    // SETUP MM METHOD
    xtp::PEwald3D3D cape = xtp::PEwald3D3D(top, xjob.getPolarTop(), _options,
		thread->getLogger());
	if (_pdb_check)
		cape.WriteDensitiesPDB(xjob.getTag()+".densities.pdb");

    // SETUP QMAPE
    QMAPEMachine machine = QMAPEMachine(&xjob, &cape, _options, "options.qmape");
    machine.setLog(thread->getLogger());
    
    // EVALUATE: ITERATE UNTIL CONVERGED
    machine.Evaluate(&xjob);

    // GENERATE OUTPUT AND FORWARD TO PROGRESS OBSERVER (RETURN)
    xtp::Job::JobResult jres = xtp::Job::JobResult();
    jres.setOutput(xjob.getInfoLine());
    jres.setStatus(xtp::Job::COMPLETE);

    return jres;
}



    
}}

#endif /* __QMAPE__H */
