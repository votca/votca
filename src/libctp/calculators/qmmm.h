#ifndef __QMMMCALC__H
#define	__QMMMCALC__H

#include <votca/ctp/parallelxjobcalc.h>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
#include <votca/ctp/xinteractor.h>
#include <votca/ctp/qmmachine.h>
#include <boost/format.hpp>


using boost::format;

namespace votca { namespace ctp {

   
class QMMM : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{

public:

    QMMM() {};
   ~QMMM() {};
   
    string          Identify() { return "QMMM"; }
    void            Initialize(Topology *, Property *);

    void            CustomizeLogger(QMThread *thread);
    void            PreProcess(Topology *top);
    Job::JobResult  EvalJob(Topology *top, Job *job, QMThread *thread);
    void            PostProcess(Topology *top);
    

private:
    
    Property                       *_options;
    
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
    bool                            _induce;
    bool                            _induce_intra_pair;

    // Multipole Interaction parameters
    string                          _method;
    bool                            _useCutoff;
    double                          _cutoff1;
    double                          _cutoff2;
    
    // QM Package options
    string                          _package;
    Property                        _qmpack_opt;

    // XJob logbook (file output)
    string                          _outFile;
    bool                            _energies2File;


};

// ========================================================================== //
//                      PARALLELCALC MEMBER FUNCTIONS                         //
// ========================================================================== //


void QMMM::Initialize(Topology *top, Property *opt) {

    _options = opt;
    
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
            _xjobfile = opt->get(key+".job_file").as<string>();
        }
        else {
            throw std::runtime_error("XJob-file not set. Abort.");
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

        if (opt->exists(key+".pdb_check")) {
            _pdb_check = opt->get(key+".pdb_check").as<string>();
        }
        else { _pdb_check = ""; }

        if (opt->exists(key+".write_chk")) {
            _write_chk_suffix = opt->get(key+".write_chk").as<string>();
            _write_chk = true;
        }
        else { _write_chk = false; }

        if (opt->exists(key+".format_chk")) {
            _chk_format = opt->get(key+".format_chk").as<string>();
        }
        else { _chk_format = "xyz"; }

        if (opt->exists(key+".split_dpl")) {
            _chk_split_dpl = (opt->get(key+".split_dpl").as<int>() == 1) ?
                         true : false;
        }
        else { _chk_split_dpl = true; }

        if (opt->exists(key+".dpl_spacing")) {
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
            _subthreads = opt->get(key+".subthreads").as< double >();
        }
        else {
            _subthreads = 1;
        }
    
    key = "options.qmmm.qmpackage";
    
        if (opt->exists(key+".package")) {
            string package_xml = opt->get(key+".package").as< string >();
            load_property_from_xml(_qmpack_opt, package_xml.c_str());
            _package = _qmpack_opt.get("package.name").as< string >();
        }
        else {
            throw runtime_error("No QM package specified.");
        }
}


void QMMM::PreProcess(Topology *top) {

//    // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
//    cout << endl << "... ... Initialize MPS-mapper: " << flush;
//    _mps_mapper.GenerateMap(_xml_file, _emp_file, top, _XJobs);
}


void QMMM::PostProcess(Topology *top) {
    
//    // WRITE OUTPUT (PRIMARILY ENERGIE SPLITTINGS)
//    FILE *out;
//    out = fopen(this->_outFile.c_str(), "w");
//    vector<XJob*> :: iterator jit;
//    for (jit = _XJobs.begin(); jit < _XJobs.end(); ++jit) {
//        (*jit)->WriteInfoLine(out);
//    }
//    fclose(out);    
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


Job::JobResult QMMM::EvalJob(Topology *top, Job *job, QMThread *thread) {
    
    // SILENT LOGGER FOR QMPACKAGE
    Logger* log = thread->getLogger();    
    Logger* qlog = new Logger();
    qlog->setReportLevel(logWARNING);
    qlog->setMultithreading(_maverick);    
    
    // GENERATE POLAR TOPOLOGY FOR JOB
    double co1 = _cutoff1;
    double co2 = _cutoff2;    
    
//    _mps_mapper.Gen_QM_MM1_MM2(top, job, co1, co2);
//    
//    LOG(logINFO,*log)
//         << job->getPolarTop()->ShellInfoStr() << flush;
//    
//    if (tools::globals::verbose)
//    job->getPolarTop()->PrintPDB(job->getTag()+"_QM0_MM1_MM2.pdb");
//
//    // INDUCTOR, QM RUNNER, QM-MM MACHINE
//    XInductor xind = XInductor(top, _options, "options.qmmm",
//        _subthreads, _maverick);    
//    xind.setLog(thread->getLogger());
//    
//    Gaussian qmpack = Gaussian(&_qmpack_opt);
//    qmpack.setLog(qlog);
//    
//    QMMachine<Gaussian> machine = QMMachine<Gaussian>(job, &xind, &qmpack, 
//        _options, "options.qmmm", _subthreads, _maverick);
//    machine.setLog(thread->getLogger());
//    
//    // EVALUATE: ITERATE UNTIL CONVERGED
//    machine.Evaluate(job);    
//    
//    // DELIVER OUTPUT & CLEAN
//    this->LockCout();
//    cout << *thread->getLogger();
//    this->UnlockCout();
//    job->setInfoLine();
//    job->getPolarTop()->~PolarTop();
    
    return Job::JobResult(job);
}



    
}}

#endif /* __QMMM__H */