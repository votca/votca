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
#ifndef __QMMMCALC__H
#define	__QMMMCALC__H

#include <votca/xtp/parallelxjobcalc.h>
#include <votca/xtp/xmapper.h>
#include <votca/xtp/xjob.h>
#include <votca/xtp/xinductor.h>
#include <votca/xtp/xinteractor.h>
// add gwbse header of excited state support
#include <votca/xtp/gwbse.h>
// --------
#include <votca/xtp/qmmachine.h>
#include <boost/format.hpp>


using boost::format;

namespace votca { namespace xtp {

   
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
    XMpsMap                        _mps_mapper;
    
    // Control over induction-state output
    string                          _pdb_check;
    bool                            _write_chk;
    string                          _write_chk_suffix;
    bool                            _chk_split_dpl;
    double                          _chk_dpl_spacing;
    string                          _chk_format;

    
    // ======================================== //
    // INDUCTION + ENERGY EVALUATION            //
    // ======================================== //

    // Induction, subthreading (-> base class)
    //bool                            _induce;
    //bool                            _induce_intra_pair;

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
    int                             _state;

    // XJob logbook (file output)
    string                          _outFile;
    bool                            _energies2File;
    
    Property                        _options;


};

// ========================================================================== //
//                      PARALLELCALC MEMBER FUNCTIONS                         //
// ========================================================================== //


void QMMM::Initialize(Property *opt) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( opt );
    _options = *opt;
    
    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;

    _maverick = (_nThreads == 1) ? true : false;
    

    string key = "options.qmmm.multipoles";

        if ( opt->exists(key) ) {
            _xml_file = opt->get(key).as< string >();
        }

    key = "options.qmmm.control";

        if ( opt->exists(key+".job_file")) {
            _jobfile = opt->get(key+".job_file").as<string>();
        }
        else {
            throw std::runtime_error("Job-file not set. Abort.");
        }

        if ( opt->exists(key+".emp_file")) {
            _emp_file   = opt->get(key+".emp_file").as<string>();
        }
        else {
            _emp_file   = opt->get(key+".emp_file").as<string>();
        }

        if ( opt->exists(key+".output") ) {
            _outFile = opt->get(key+".output").as< string >();
            _energies2File = true;
        }
        else { _energies2File = false; }

        if ( opt->exists(key+".pdb_check")) {
            _pdb_check = opt->get(key+".pdb_check").as<string>();
        }
        else { _pdb_check = ""; }

        if ( opt->exists(key+".write_chk")) {
            _write_chk_suffix = opt->get(key+".write_chk").as<string>();
            _write_chk = true;
        }
        else { _write_chk = false; }

        if ( opt->exists(key+".format_chk")) {
            _chk_format = opt->get(key+".format_chk").as<string>();
        }
        else { _chk_format = "xyz"; }

        if ( opt->exists(key+".split_dpl")) {
            _chk_split_dpl = ( opt->get(key+".split_dpl").as<int>() == 1) ?
                         true : false;
        }
        else { _chk_split_dpl = true; }

        if ( opt->exists(key+".dpl_spacing")) {
            _chk_dpl_spacing = opt->get(key+".dpl_spacing").as<double>();
        }
        else {
            _chk_dpl_spacing = 1.0e-6;
        }


    key = "options.qmmm.coulombmethod";
    
        if ( opt->exists(key+".method") ) {
            _method = opt->get(key+".method").as< string >();
            if (_method != "cut-off" && _method != "cutoff") {
                throw runtime_error("Method " + _method + " not recognised.");
            }
        }
        else {
            _method = "cut-off";
        }
        if ( opt->exists(key+".cutoff1") ) {
            _cutoff1 = opt->get(key+".cutoff1").as< double >();
            if (_cutoff1) { _useCutoff = true; }
        }
        if ( opt->exists(key+".cutoff2") ) {
            _cutoff2 = opt->get(key+".cutoff2").as< double >();
        }
        else {
            _cutoff2 = _cutoff1;
        }
        if ( opt->exists(key+".subthreads") ) {
            _subthreads = opt->get(key+".subthreads").as< int >();
        }
        else {
            _subthreads = 1;
        }
    
    key = "options.qmmm.qmpackage";
    
        if ( opt->exists(key+".package")) {
            string package_xml = opt->get(key+".package").as< string >();
            load_property_from_xml(_qmpack_opt, package_xml.c_str());
            _package = _qmpack_opt.get("package.name").as< string >();
        }
        else {
            throw runtime_error("No QM package specified.");
        }
    
    
    // GWBSE options, depending on whether it is there, decide for ground
    // or excited state QM/MM
    key = "options.qmmm.gwbse";
    
    if ( opt->exists(key)) { 
        cout << " Excited state QM/MM " << endl;
        
         if ( opt->exists(key+".gwbse_options")) {
            string gwbse_xml = opt->get(key+".gwbse_options").as< string >();
            load_property_from_xml(_gwbse_opt, gwbse_xml.c_str());
            // _gwbse = _gwbse_opt.get("package.name").as< string >();
        }
        else {
            throw runtime_error("GWBSE options not specified.");
        }
        
        _state = opt->get(key+".state").as< int >();
        
        
    } else {
        cout << " Ground state QM/MM " << endl;
    }
    

    
    
    //cout << TXT << _options;
    
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
            LOG(logERROR,*(thread->getLogger()))
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
    Logger* qlog = new Logger();
    qlog->setReportLevel(logDEBUG);
    qlog->setMultithreading(_maverick);
    qlog->setPreface(logINFO,    (format("\nQ%1$02d ... ...") % thread->getId()).str());
    qlog->setPreface(logERROR,   (format("\nQ%1$02d ERR ...") % thread->getId()).str());
    qlog->setPreface(logWARNING, (format("\nQ%1$02d WAR ...") % thread->getId()).str());
    qlog->setPreface(logDEBUG,   (format("\nQ%1$02d DBG ...") % thread->getId()).str());      
    
    // CREATE XJOB FROM JOB INPUT STRING
    LOG(logINFO,*log)
        << "Job input = " << job->getInput().as<string>() << flush;
    XJob xjob = this->ProcessInputString(job, top, thread);  
    
    // GENERATE POLAR TOPOLOGY FOR JOB
    double co1 = _cutoff1;
    double co2 = _cutoff2;    
    _mps_mapper.Gen_QM_MM1_MM2(top, &xjob, co1, co2, thread);
    
    LOG(logINFO,*log)
         << xjob.getPolarTop()->ShellInfoStr() << flush;
    
    if (tools::globals::verbose)
    xjob.getPolarTop()->PrintPDB(xjob.getTag()+"_QM0_MM1_MM2.pdb");

    // INDUCTOR, QM RUNNER, QM-MM MACHINE
    XInductor xind = XInductor(top, &_options, "options.qmmm",
        _subthreads, _maverick);    
    xind.setLog(thread->getLogger());
    
    //Gaussian qmpack = Gaussian(&_qmpack_opt);

    // get the corresponding object from the QMPackageFactory
    QMPackage *qmpack =  QMPackages().Create( _package );
    
    qmpack->Initialize( &_qmpack_opt );
    
    qmpack->setLog(qlog);
    
    QMMachine<QMPackage> machine = QMMachine<QMPackage>(&xjob, &xind, qmpack, 
        &_options, "options.qmmm", _subthreads, _maverick);
    machine.setLog(thread->getLogger());
    
    // EVALUATE: ITERATE UNTIL CONVERGED
    machine.Evaluate(&xjob);    
    
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
        LOG(logERROR,*log) << xind.getError() << flush;
    }
    
    return jres;
}



    
}}

#endif /* __QMMM__H */
