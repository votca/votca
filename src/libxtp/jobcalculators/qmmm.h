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

#include <votca/ctp/parallelxjobcalc.h>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
#include <votca/ctp/xinteractor.h>
// add gwbse header of excited state support
#include <votca/xtp/gwbse.h>
// --------
#include <votca/xtp/qmmachine.h>
#include <boost/format.hpp>


using boost::format;

namespace votca { namespace xtp {


class QMMM : public ctp::ParallelXJobCalc< vector<ctp::Job*>, ctp::Job*, ctp::Job::JobResult >
{

public:

    QMMM() {};
   ~QMMM() {};

    string          Identify() { return "qmmm"; }
    void            Initialize(Property *);

    void            CustomizeLogger(ctp::QMThread *thread);
    void            PreProcess(ctp::Topology *top);
    ctp::Job::JobResult  EvalJob(ctp::Topology *top, ctp::Job *job, ctp::QMThread *thread);
    ctp::XJob            ProcessInputString(ctp::Job *job, ctp::Topology *top, ctp::QMThread *thread);


private:


    // ======================================== //
    // MULTIPOLE ALLOCATION, XJOBS, ADD. OUTPUT //
    // ======================================== //

    // Polar-site mapping
    string                         _emp_file;
    string                         _xml_file;
    ctp::XMpsMap                        _mps_mapper;

    // Control over induction-state output
    string                          _pdb_check;
    bool                            _write_chk;
    string                          _write_chk_suffix;
    string                          _chk_format;


    // ======================================== //
    // INDUCTION + ENERGY EVALUATION            //
    // ======================================== //

  
    // Multipole Interaction parameters
    string                          _method;
    bool                            _useCutoff;
    double                          _cutoff1;
    double                          _cutoff2;

    // QM Package options
    string                          _package;
    Property                        _qmpack_opt;

    // GWBSE options
    string                          _gwbse;
    Property                        _gwbse_opt;

    // XJob logbook (file output)
    string                          _outFile;
    bool                            _energies2File;

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


    string key = "options."+Identify()+".mapping";

        if ( options->exists(key) ) {
            _xml_file = options->get(key).as< string >();
        }

    key = "options."+Identify();

        if ( options->exists(key+".job_file")) {
            _jobfile = options->get(key+".job_file").as<string>();
        }
        else {
            throw std::runtime_error("Job-file not set. Abort.");
        }

        if ( options->exists(key+".emp_file")) {
            _emp_file   = options->get(key+".emp_file").as<string>();
        }
        else {
            _emp_file   = options->get(key+".emp_file").as<string>();
        }

        if ( options->exists(key+".output") ) {
            _outFile = options->get(key+".output").as< string >();
            _energies2File = true;
        }
        else { _energies2File = false; }

        if ( options->exists(key+".pdb_check")) {
            _pdb_check = options->get(key+".pdb_check").as<string>();
        }
        else { _pdb_check = ""; }

        if ( options->exists(key+".write_chk")) {
            _write_chk_suffix = options->get(key+".write_chk").as<string>();
            _write_chk = true;
        }
        else { _write_chk = false; }

        if ( options->exists(key+".format_chk")) {
            _chk_format = options->get(key+".format_chk").as<string>();
        }
        else { _chk_format = "xyz"; }

        


    key = "options."+Identify()+".coulombmethod";

        if ( options->exists(key+".method") ) {
            _method = options->get(key+".method").as< string >();
            if (_method != "cut-off" && _method != "cutoff") {
                throw runtime_error("Method " + _method + " not recognised.");
            }
        }
        else {
            _method = "cut-off";
        }
        if ( options->exists(key+".cutoff1") ) {
            _cutoff1 = options->get(key+".cutoff1").as< double >();
            if (_cutoff1) { _useCutoff = true; }
        }
        if ( options->exists(key+".cutoff2") ) {
            _cutoff2 = options->get(key+".cutoff2").as< double >();
        }
        else {
            _cutoff2 = _cutoff1;
        }
    if(_cutoff1>_cutoff2){
        throw runtime_error("Cutoff1 must be smaller or equal Cutoff2");
    }
        if ( options->exists(key+".subthreads") ) {
            _subthreads = options->get(key+".subthreads").as< int >();
        }
        else {
            _subthreads = 1;
        }

    key = "options."+Identify();

        if ( options->exists(key+".dftpackage")) {
            string package_xml = options->get(key+".dftpackage").as< string >();
            load_property_from_xml(_qmpack_opt, package_xml.c_str());
            _package = _qmpack_opt.get("package.name").as< string >();
        }
        else {
            throw runtime_error("No QM package specified.");
        }


    // GWBSE options, depending on whether it is there, decide for ground
    // or excited state QM/MM
    key = "options."+Identify()+".gwbse";

    if ( options->exists(key)) {
         if ( options->exists(key+".gwbse_options")) {
            string gwbse_xml = options->get(key+".gwbse_options").as< string >();
            load_property_from_xml(_gwbse_opt, gwbse_xml.c_str());
        }
        else {
            throw runtime_error("GWBSE options not specified.");
        }

    }

    // register all QM packages (Gaussian, turbomole, etc))
    QMPackageFactory::RegisterAll();

}


void QMMM::PreProcess(ctp::Topology *top) {

    // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
    cout << endl << "... ... Initialize MPS-mapper: " << flush;
    _mps_mapper.GenerateMap(_xml_file, _emp_file, top);
}


void QMMM::CustomizeLogger(ctp::QMThread *thread) {

    // CONFIGURE LOGGER
    ctp::Logger* log = thread->getLogger();
    log->setReportLevel(ctp::logDEBUG);
    log->setMultithreading(_maverick);

    log->setPreface(ctp::logINFO,    (format("\nT%1$02d ... ...") % thread->getId()).str());
    log->setPreface(ctp::logERROR,   (format("\nT%1$02d ERR ...") % thread->getId()).str());
    log->setPreface(ctp::logWARNING, (format("\nT%1$02d WAR ...") % thread->getId()).str());
    log->setPreface(ctp::logDEBUG,   (format("\nT%1$02d DBG ...") % thread->getId()).str());
}


// ========================================================================== //
//                            QMMM MEMBER FUNCTIONS                           //
// ========================================================================== //


ctp::XJob QMMM::ProcessInputString(ctp::Job *job, ctp::Topology *top, ctp::QMThread *thread) {

    string input = job->getInput().as<string>();
    vector<ctp::Segment*> qmSegs;
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

        ctp::Segment *seg = top->getSegment(segId);
        if (seg->getName() != segName) {
            CTP_LOG(ctp::logERROR,*(thread->getLogger()))
                << "ERROR: Seg " << segId << ":" << seg->getName() << " "
                << " maltagged as " << segName << ". Skip job ..." << flush;
            throw std::runtime_error("Input does not match topology.");
        }

        qmSegs.push_back(seg);
        qmSegMps.push_back(mpsFile);
    }

    return ctp::XJob(job->getId(), job->getTag(), qmSegs, qmSegMps, top);
}


ctp::Job::JobResult QMMM::EvalJob(ctp::Topology *top, ctp::Job *job, ctp::QMThread *thread) {

    // SILENT LOGGER FOR QMPACKAGE
    ctp::Logger* log = thread->getLogger();
    ctp::Logger* qlog = new ctp::Logger();
    qlog->setReportLevel(ctp::logDEBUG);
    qlog->setMultithreading(_maverick);
    qlog->setPreface(ctp::logINFO,    (format("\nQ%1$02d ... ...") % thread->getId()).str());
    qlog->setPreface(ctp::logERROR,   (format("\nQ%1$02d ERR ...") % thread->getId()).str());
    qlog->setPreface(ctp::logWARNING, (format("\nQ%1$02d WAR ...") % thread->getId()).str());
    qlog->setPreface(ctp::logDEBUG,   (format("\nQ%1$02d DBG ...") % thread->getId()).str());

    // CREATE XJOB FROM JOB INPUT STRING
    CTP_LOG(ctp::logINFO,*log)
        << "Job input = " << job->getInput().as<string>() << flush;
    ctp::XJob xjob = this->ProcessInputString(job, top, thread);

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


    CTP_LOG(ctp::logINFO,*log)
         << xjob.getPolarTop()->ShellInfoStr() << flush;

    if (tools::globals::verbose){
        xjob.getPolarTop()->PrintPDB(xjob.getTag()+"_QM0_MM1_MM2.pdb");
    }
    // INDUCTOR, QM RUNNER, QM-MM MACHINE
    ctp::XInductor xind = ctp::XInductor(top, &_options, "options."+Identify(),
        _subthreads, _maverick);
    xind.setLog(thread->getLogger());

    //Gaussian qmpack = Gaussian(&_qmpack_opt);

    // get the corresponding object from the QMPackageFactory
    QMPackage *qmpack =  QMPackages().Create( _package );

    qmpack->Initialize( _qmpack_opt );

    qmpack->setLog(qlog);
    

    QMMachine<QMPackage> machine = QMMachine<QMPackage>(&xjob, &xind, qmpack,
        &_options, "options."+Identify(), _subthreads, _maverick);
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
    ctp::Job::JobResult jres = ctp::Job::JobResult();
    jres.setOutput(xjob.getInfoLine());
    jres.setStatus(ctp::Job::COMPLETE);

    if (!xind.hasConverged()) {
        jres.setStatus(ctp::Job::FAILED);
        jres.setError(xind.getError());
        CTP_LOG(ctp::logERROR,*log) << xind.getError() << flush;
    }
    if(error!=0){
        jres.setStatus(ctp::Job::FAILED);
    }

    return jres;
}




}}

#endif /* __QMMM__H */
