#ifndef XQMULTIPOLE_H
#define XQMULTIPOLE_H


#include <votca/xtp/parallelxjobcalc.h>
#include <votca/xtp/xmapper.h>
#include <votca/xtp/xjob.h>
#include <votca/xtp/xinductor.h>
#include <votca/xtp/xinteractor.h>
#include <votca/xtp/logger.h>
#include <boost/format.hpp>


using boost::format;


namespace votca { namespace xtp {

    
class XQMP : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{

public:

    XQMP() {};
   ~XQMP() {};
   
    string          Identify() { return "xqmultipole"; }
    void            Initialize(Property *);

    void            CustomizeLogger(QMThread *thread);
    void            PreProcess(Topology *top);
    Job::JobResult  EvalJob(Topology *top, Job *job, QMThread *thread);
    void            PostProcess(Topology *top) { ; }
    
    XJob            ProcessInputString(Job *, Topology *, QMThread *);
    

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
    bool                            _pdb_check;
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



};

// ========================================================================== //
//                        XMULTIPOLE MEMBER FUNCTIONS                         //
// ========================================================================== //


void XQMP::Initialize(Property *opt) {

    _options = opt;
    
    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;

    _maverick = (_nThreads == 1) ? true : false;
    

    string key = "options.xqmultipole.multipoles";

        if ( opt->exists(key) ) {
            _xml_file = opt->get(key).as< string >();
        }

    key = "options.xqmultipole.control";

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

        if (opt->exists(key+".pdb_check")) {
            _pdb_check = opt->get(key+".pdb_check").as<bool>();
        }
        else { _pdb_check = false; }

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


    key = "options.xqmultipole.coulombmethod";
    
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
}


void XQMP::PreProcess(Topology *top) {

    // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
    cout << endl << "... ... Initialize MPS-mapper: " << flush;
    _mps_mapper.GenerateMap(_xml_file, _emp_file, top);
}


void XQMP::CustomizeLogger(QMThread *thread) {
    
    // CONFIGURE LOGGER
    Logger* log = thread->getLogger();
    log->setReportLevel(logDEBUG);
    log->setMultithreading(_maverick);

    log->setPreface(logINFO,    (format("\nT%1$02d INF ...") % thread->getId()).str());
    log->setPreface(logERROR,   (format("\nT%1$02d ERR ...") % thread->getId()).str());
    log->setPreface(logWARNING, (format("\nT%1$02d WAR ...") % thread->getId()).str());
    log->setPreface(logDEBUG,   (format("\nT%1$02d DBG ...") % thread->getId()).str());        
}



// ========================================================================== //
//                           JOBXQMP MEMBER FUNCTIONS                         //
// ========================================================================== //


XJob XQMP::ProcessInputString(Job *job, Topology *top, QMThread *thread) {
    
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


Job::JobResult XQMP::EvalJob(Topology *top, Job *job, QMThread *thread) {
    
    Logger *log = thread->getLogger();
    
    LOG(logINFO,*log)
        << "Job input = " << job->getInput().as<string>() << flush;
    
    // CREATE XJOB FROM JOB INPUT STRING
    XJob xjob = this->ProcessInputString(job, top, thread);
    
    
    // GENERATE POLAR TOPOLOGY
    double co1 = _cutoff1;
    double co2 = _cutoff2;    
    
    _mps_mapper.Gen_QM_MM1_MM2(top, &xjob, co1, co2, thread);
    
    LOG(logINFO,*log)
         << xjob.getPolarTop()->ShellInfoStr() << flush;
    
    if (tools::globals::verbose || _pdb_check)
    xjob.getPolarTop()->PrintPDB(xjob.getTag()+"_QM0_MM1_MM2.pdb");

    // CALL MAGIC INDUCTOR         
    XInductor inductor = XInductor(top, _options, "options.xqmultipole",
                                   _subthreads, _maverick);
    inductor.setLog(thread->getLogger());
    inductor.Evaluate(&xjob);
    

    // SAVE INDUCTION STATE
    if (_write_chk) {
        
        string format    = _chk_format;
        string dotsuffix = (format == "gaussian") ? ".com" : ".xyz";
        string outstr    = xjob.getTag()+_write_chk_suffix+dotsuffix;
        
        bool split       = _chk_split_dpl;
        double space     = _chk_dpl_spacing;
        
        xjob.getPolarTop()->PrintInduState(outstr, format, split, space);        
    }

    // JOT INFO STRING & CLEAN POLAR TOPOLOGY
    xjob.setInfoLine(true,false);
    
    
    // GENERATE OUTPUT AND FORWARD TO PROGRESS OBSERVER (RETURN)
    Property output = xjob.GenerateOutputProperty();
    Job::JobResult jres = Job::JobResult();
    jres.setOutput(output);
    jres.setStatus(Job::COMPLETE);
    
    if (!inductor.hasConverged()) {
        jres.setStatus(Job::FAILED);
        jres.setError(inductor.getError());
        LOG(logERROR,*log) << inductor.getError() << flush;
    }
    
    return jres;
}


}}

#endif

