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
#ifndef __QMMMCALC__H
#define	__QMMMCALC__H

#include <votca/xtp/parallelxjobcalc.h>
#include <votca/xtp/xmapper.h>
#include <votca/xtp/xjob.h>
#include <votca/xtp/xinductor.h>
#include <votca/xtp/xinteractor.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmmachine.h>
#include <votca/xtp/polartop.h>
#include <boost/format.hpp>

namespace votca { namespace xtp {
using boost::format;

class QMMM : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{

public:

    QMMM() {};
   ~QMMM() {};

    string          Identify() { return "qmmm"; }
    void            Initialize(Property *);

    void            CustomizeLogger(QMThread *thread);
    void            PreProcess(Topology *top);
    Job::JobResult  EvalJob(Topology *top, Job *job, QMThread *thread);
    XJob            ProcessInputString(Job *job, Topology *top, QMThread *thread);


private:


    // ======================================== //
    // MULTIPOLE ALLOCATION, XJOBS, ADD. OUTPUT //
    // ======================================== //

    // Polar-site mapping
    string                         _emp_file;
    string                         _xml_file;
    XMpsMap                   _mps_mapper;

    // ======================================== //
    // INDUCTION + ENERGY EVALUATION            //
    // ======================================== //
  
    // Multipole Interaction parameters
    string                          _method;
    double                          _cutoff1;
    double                          _cutoff2;

    // QM Package options
    string                          _package;
    Property                        _qmpack_opt;

    // GWBSE options
    string                          _gwbse;
    Property                        _gwbse_opt;

    Property                        _options;


};

// ========================================================================== //
//                      PARALLELCALC MEMBER FUNCTIONS                         //
// ========================================================================== //


void QMMM::Initialize(Property *options) {

    // update options with the VOTCASHARE defaults
    UpdateWithDefaults( options, "xtp" );
    _options = *options;

    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;

    _maverick = (_nThreads == 1) ? true : false;


    string key = "options."+Identify();
    _xml_file= options->ifExistsReturnElseThrowRuntimeError<string>(key+".mapping"); 
    _jobfile = options->ifExistsReturnElseThrowRuntimeError<string>(key+".job_file"); 
    _emp_file   = options->ifExistsReturnElseThrowRuntimeError<string>(key+".emp_file"); 
    
    key = "options."+Identify()+".coulombmethod";

    std::vector<string> choices={"cutoff","cut-off"};
    _method = options->ifExistsAndinListReturnElseThrowRuntimeError<string>(key+".method",choices);
    _cutoff1 = options->ifExistsReturnElseThrowRuntimeError<double>(key+".cutoff1");
    _cutoff2 = options->ifExistsReturnElseReturnDefault<double>(key+".cutoff2",_cutoff1);
    
    if(_cutoff1>_cutoff2){
        throw runtime_error("Cutoff1 must be smaller or equal Cutoff2");
    }
 
    _subthreads = options->ifExistsReturnElseReturnDefault<int>(key+".subthreads",1);
      

    key = "options."+Identify();

     string package_xml = options->ifExistsReturnElseThrowRuntimeError<string>(key+".dftpackage");
    load_property_from_xml(_qmpack_opt, package_xml.c_str());
    _package = _qmpack_opt.get("package.name").as< string >();
       


    // GWBSE options, depending on whether it is there, decide for ground
    // or excited state QM/MM
    key = "options."+Identify()+".gwbse";

    if ( options->exists(key)) {
            string gwbse_xml = options->ifExistsReturnElseThrowRuntimeError<string>(key+".gwbse_options");
            load_property_from_xml(_gwbse_opt, gwbse_xml.c_str());
    }

    // register all QM packages (Gaussian, turbomole, etc))
    QMPackageFactory::RegisterAll();

}


void QMMM::PreProcess(Topology *top) {

    // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
    cout << endl << "... ... Initialize MPS-mapper: " << flush;
    _mps_mapper.GenerateMap(_xml_file, _emp_file, top);
}


void QMMM::CustomizeLogger(QMThread *thread) {

    // CONFIGURE LOGGER
    Logger* log = thread->getLogger();
    log->setReportLevel(logDEBUG);
    log->setMultithreading(_maverick);

    log->setPreface(logINFO,    (format("\nT%1$02d ... ...") % thread->getId()).str());
    log->setPreface(logERROR,   (format("\nT%1$02d ERR ...") % thread->getId()).str());
    log->setPreface(logWARNING, (format("\nT%1$02d WAR ...") % thread->getId()).str());
    log->setPreface(logDEBUG,   (format("\nT%1$02d DBG ...") % thread->getId()).str());
}


// ========================================================================== //
//                            QMMM MEMBER FUNCTIONS                           //
// ========================================================================== //


XJob QMMM::ProcessInputString(Job *job, Topology *top, QMThread *thread) {

    string input = job->getInput().as<string>();
    vector<Segment*> qmSegs;
    vector<string>   qmSegMps;
    vector<string> split;
    Tokenizer toker(input, " \t\n");
    toker.ToVector(split);

    for (unsigned int i = 0; i < split.size(); ++i) {

        string id_seg_mps = split[i];
        vector<string> split_id_seg_mps;
        Tokenizer toker(id_seg_mps, ":");
        toker.ToVector(split_id_seg_mps);

        int segId = boost::lexical_cast<int>(split_id_seg_mps[0]);
        string segName = split_id_seg_mps[1];
        string mpsFile = split_id_seg_mps[2];

        Segment *seg = top->getSegment(segId);
        if (seg->getName() != segName) {
            XTP_LOG(logERROR,*(thread->getLogger()))
                << "ERROR: Seg " << segId << ":" << seg->getName() << " "
                << " maltagged as " << segName << ". Skip job ..." << flush;
            throw std::runtime_error("Input does not match topology.");
        }

        qmSegs.push_back(seg);
        qmSegMps.push_back(mpsFile);
    }

    return XJob(job->getId(), job->getTag(), qmSegs, qmSegMps, top);
}


Job::JobResult QMMM::EvalJob(Topology *top, Job *job, QMThread *thread) {

    // SILENT LOGGER FOR QMPACKAGE
    Logger* log = thread->getLogger();
    Logger qlog;
    qlog.setReportLevel(logDEBUG);
    qlog.setMultithreading(_maverick);
    qlog.setPreface(logINFO,    (format("\nQ%1$02d ... ...") % thread->getId()).str());
    qlog.setPreface(logERROR,   (format("\nQ%1$02d ERR ...") % thread->getId()).str());
    qlog.setPreface(logWARNING, (format("\nQ%1$02d WAR ...") % thread->getId()).str());
    qlog.setPreface(logDEBUG,   (format("\nQ%1$02d DBG ...") % thread->getId()).str());

    // CREATE XJOB FROM JOB INPUT STRING
    XTP_LOG(logINFO,*log)
        << "Job input = " << job->getInput().as<string>() << flush;
    XJob xjob = this->ProcessInputString(job, top, thread);

    // GENERATE POLAR TOPOLOGY FOR JOB
    double co1 = _cutoff1;
    double co2 = _cutoff2;
    _mps_mapper.Gen_QM_MM1_MM2(top, &xjob, co1, co2, thread);

    const matrix box=xjob.getTop()->getBox();
    //check if box is non orthogonal

    double min=box.get(0,0);
    if(min>box.get(1,1)){min=box.get(1,1);}
    if(min>box.get(2,2)){min=box.get(2,2);}

    if(_cutoff2>0.5*min){
        throw runtime_error((format("Cutoff is larger than half the box size. Maximum allowed cutoff is %1$1.1f - molecule extension.") % (0.5*min)).str());
    }


    XTP_LOG(logINFO,*log)
         << xjob.getPolarTop()->ShellInfoStr() << flush;

    if (tools::globals::verbose){
        xjob.getPolarTop()->PrintPDB(xjob.getTag()+"_QM0_MM1_MM2.pdb");
    }
    // INDUCTOR, QM RUNNER, QM-MM MACHINE
    XInductor xind = XInductor(top, &_options, "options."+Identify(),
        _subthreads, _maverick);
    xind.setLog(thread->getLogger());

    // get the corresponding object from the QMPackageFactory
    QMPackage *qmpack =  QMPackages().Create( _package );
    qmpack->Initialize( _qmpack_opt );
    qmpack->setLog(&qlog);
    

    QMMachine machine = QMMachine(&xjob, &xind, qmpack,
        &_options, "options."+Identify());
    machine.setLog(thread->getLogger());

    // EVALUATE: ITERATE UNTIL CONVERGED
    int error=machine.Evaluate(&xjob);

    // DESTROY QMPackage
    delete qmpack;

    // DELIVER OUTPUT & CLEAN
    this->LockCout();
    cout << *thread->getLogger();
    this->UnlockCout();

    // JOT INFO STRING & CLEAN POLAR TOPOLOGY
    xjob.setInfoLine(true,true);

    // GENERATE OUTPUT AND FORWARD TO PROGRESS OBSERVER (RETURN)
    Job::JobResult jres = Job::JobResult();
    jres.setOutput(xjob.getInfoLine());
    jres.setStatus(Job::COMPLETE);

    if (!xind.hasConverged()) {
        jres.setStatus(Job::FAILED);
        jres.setError(xind.getError());
        XTP_LOG(logERROR,*log) << xind.getError() << flush;
    }
    if(error!=0){
        jres.setStatus(Job::FAILED);
    }

    return jres;
}




}}

#endif /* __QMMM__H */
