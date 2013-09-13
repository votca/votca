#ifndef VOTCA_CTP_EWALD_H
#define VOTCA_CTP_EWALD_H


#include <votca/ctp/parallelxjobcalc.h>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
#include <votca/ctp/xinteractor.h>
#include <votca/ctp/ewald2d.h>
#include <votca/ctp/ewald3d.h>
#include <votca/ctp/pewald3d.h>
#include <votca/ctp/logger.h>
#include <boost/format.hpp>


using boost::format;


namespace votca { namespace ctp {


template<class EwaldMethod>
class Ewald : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{

public:

    Ewald() {};
   ~Ewald() {};
   
    string          Identify() { return "ewald"; }
    void            Initialize(Topology *, Property *);

    void            PreProcess(Topology *top);
    Job::JobResult  EvalJob(Topology *top, Job *job, QMThread *thread);
    void            PostProcess(Topology *top) { ; }
    
    XJob            ProcessInputString(const Job *, Topology *, QMThread *);
    

private:
    
    Property                       *_options;
    
    // ======================================== //
    // MULTIPOLE ALLOCATION, XJOBS, ADD. OUTPUT //
    // ======================================== //

    // Polar-site mapping
    string                         _mps_table;
    string                         _xml_file;
    XMpsMap                        _mps_mapper;
    bool                           _pdb_check;

    // Multipole Interaction parameters
    string                          _method;
    double                          _cutoff;



};

// ========================================================================== //
//                           EWALD MEMBER FUNCTIONS                           //
// ========================================================================== //


template<class EwaldMethod>
void Ewald<EwaldMethod>::Initialize(Topology *top, Property *opt) {

    _options = opt;
    
    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;

    _maverick = (_nThreads == 1) ? true : false;
    

    string key = "options.ewald.multipoles";

        if ( opt->exists(key) ) {
            _xml_file = opt->get(key).as< string >();
        }

    key = "options.ewald.control";

        if ( opt->exists(key+".job_file")) {
            _jobfile = opt->get(key+".job_file").as<string>();
        }
        else {
            throw std::runtime_error("Job-file not set. Abort.");
        }

        if ( opt->exists(key+".mps_table")) {
            _mps_table = opt->get(key+".mps_table").as<string>();
        }
        else {
            _mps_table = opt->get(key+".emp_file").as<string>();
        }

        if (opt->exists(key+".pdb_check")) {
            _pdb_check = opt->get(key+".pdb_check").as<bool>();
        }
        else { _pdb_check = false; }


    key = "options.ewald.coulombmethod";
    
        if ( opt->exists(key+".method") ) {
            _method = opt->get(key+".method").as< string >();
            if (_method != "ewald" && _method != "Ewald") {
                throw runtime_error("Method " + _method + " not recognised.");
            }
        }
        else {
            _method = "ewald";
        }
        if ( opt->exists(key+".cutoff") ) {
            _cutoff = opt->get(key+".cutoff").as< double >();
        }
        else {
            throw runtime_error("No real-space cut-off provided.");
        }
    
    _mps_mapper.setEstaticsOnly(true);
    return;
}


template<class EwaldMethod>
void Ewald<EwaldMethod>::PreProcess(Topology *top) {

    // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
    cout << endl << "... ... Initialize MPS-mapper: " << flush;
    _mps_mapper.GenerateMap(_xml_file, _mps_table, top);
    return;
}


template<class EwaldMethod>
XJob Ewald<EwaldMethod>::ProcessInputString(const Job *job, Topology *top, QMThread *thread) {
    
    string input = job->getInput();
    vector<Segment*> qmSegs;
    vector<string>   qmSegMps;
    vector<string> split;
    Tokenizer toker(input, " \t\n");
    toker.ToVector(split);

    for (int i = 0; i < split.size(); ++i) {
                
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


template<class EwaldMethod>
Job::JobResult Ewald<EwaldMethod>::EvalJob(Topology *top, Job *job, QMThread *thread) {
    
    Logger *log = thread->getLogger();    
    LOG(logINFO,*log)
        << "Job input = " << job->getInput() << flush;
    
    // CREATE XJOB FROM JOB INPUT STRING
    XJob xjob = this->ProcessInputString(job, top, thread);    
    
    // GENERATE POLAR TOPOLOGY
    _mps_mapper.Gen_FGC_FGN_BGN(top, &xjob, thread);
    
    // CALL EWALD MAGIC
    EwaldMethod ewaldnd = EwaldMethod(top, xjob.getPolarTop(), _options, thread->getLogger());
    if (tools::globals::verbose || _pdb_check)
        ewaldnd.WriteDensitiesPDB(xjob.getTag()+"_ew_densities.pdb");
    ewaldnd.Evaluate();
    
    // GENERATE OUTPUT AND FORWARD TO PROGRESS OBSERVER (RETURN)
    Property output = ewaldnd.GenerateOutputString();
    Job::JobResult jres = Job::JobResult();
    jres.setOutput(output);
    jres.setStatus(Job::COMPLETE);
    
    if (!ewaldnd.Converged()) {
        jres.setStatus(Job::FAILED);
        jres.setError(ewaldnd.GenerateErrorString());
        LOG(logERROR,*log) << ewaldnd.GenerateErrorString() << flush;
    }    
    
    return jres;
}


}}

#endif

