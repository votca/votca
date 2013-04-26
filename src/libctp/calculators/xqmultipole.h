#ifndef XQMULTIPOLE_H
#define XQMULTIPOLE_H


#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
#include <votca/ctp/xinteractor.h>


// TODO JobXQMP requires destructor for privatized polar sites.

namespace votca { namespace ctp {

    
class XQMP : public QMCalculator
{

public:

    XQMP() {};
   ~XQMP() {};

    string          Identify() { return "XQMultipole"; }
    void            Initialize(Topology *, Property *);
    bool            EvaluateFrame(Topology *top);

    XJob           *RequestNextJob(int id, Topology *top);
    void            LockCout() { _coutMutex.Lock(); }
    void            UnlockCout() { _coutMutex.Unlock(); }
    
    // ======================================== //
    // XJOB OPERATOR (THREAD)                   //
    // ======================================== //

    class JobXQMP : public Thread
    {
    public:

        JobXQMP(int id,   Topology *top, XQMP *master)
              : _id(id), _top(top),     _master(master) {};
       ~JobXQMP() {};

        int         getId()             { return _id; }
        void        setId(int id)       { _id = id; }
        void        InitSlotData(Topology *top) { ; }
        void        Run(void);
        void        EvalJob(Topology *top, XJob *job);

    public:

        int         _id;
        Topology   *_top;
        XQMP       *_master;
        XJob       *_job;
        XInductor   _inductor;

    };

private:

    // ======================================== //
    // MULTIPOLE ALLOCATION, XJOBS, ADD. OUTPUT //
    // ======================================== //

    // Polar-site mapping
    string                         _job_file;
    string                         _emp_file;
    string                         _xml_file;
    XMpsMap                        _mps_mapper;
    
    // XJobs and thread management
    vector<XJob*>                  _XJobs;
    vector<XJob*>::iterator        _nextXJob;
    Mutex                          _nextJobMutex;
    Mutex                          _coutMutex;
    Mutex                          _logMutex;
    bool                           _maverick;
    
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

    // Induction switches, subthreading
    bool                            _induce;
    bool                            _induce_intra_pair;
    int                             _subthreads;

    // Multipole Interaction parameters
    bool                            _useCutoff;
    double                          _cutoff1;
    double                          _cutoff2;
    bool                            _useExp;
    double                          _aDamp;
    bool                            _useScaling;
    vector<double>                  _scale1;

    // Convergence parameters for induction
    float                           _wSOR_N;
    float                           _wSOR_C;
    double                          _epsTol;
    int                             _maxIter;

    // XJob logbook (file output)
    string                          _outFile;
    bool                            _energies2File;


};

// ========================================================================== //
//                        XMULTIPOLE MEMBER FUNCTIONS                         //
// ========================================================================== //


void XQMP::Initialize(Topology *top, Property *opt) {

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
            _job_file   = opt->get(key+".job_file").as<string>();
        }
        else {
            _job_file   = opt->get(key+".job_file").as<string>();
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


    key = "options.xqmultipole.tholemodel";

        if ( opt->exists(key+".induce") ) {
            int induce = opt->get(key+".induce").as< int >();
            _induce = (induce == 0) ? false : true;
        }
        else { _induce = true; }

        if ( opt->exists(key+".induce_intra_pair") ) {
            int induce = opt->get(key+".induce_intra_pair").as< int >();
            _induce_intra_pair = (induce == 0) ? false : true;
        }
        else { _induce_intra_pair = true; }

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
        if ( opt->exists(key+".exp_damp") ) {
            _aDamp = opt->get(key+".exp_damp").as< double >();
            if (_aDamp) { _useExp = true; }
        }
        else {
            cout << endl << "... ... WARNING: No sharpness parameter supplied";
            cout << endl << "... ... ... Using default a = 0.39";
            _aDamp = 0.39;
            _useExp = true;
        }
        if ( opt->exists(key+".scaling") ) {
            _scale1 = opt->get(key+".scaling").as< vector<double> >();
            if (0 < _scale1.size() && _scale1.size() < 4) {
                _useScaling = true; }
            else {
                _useScaling = false;
                cout << endl << "... ... WARNING: 1-N SCALING SWITCHED OFF"; }
        }

    key = "options.xqmultipole.convergence";

        if ( opt->exists(key+".wSOR_N") ) {
            _wSOR_N = opt->get(key+".wSOR_N").as< float >();
        }
        else { _wSOR_N = 0.75; }
        if ( opt->exists(key+".wSOR_C") ) {
            _wSOR_C = opt->get(key+".wSOR_C").as< float >();
        }
        else { _wSOR_C = 0.75; }

        if ( opt->exists(key+".max_iter") ) {
            _maxIter = opt->get(key+".max_iter").as< int >();
        }
        else { _maxIter = 512; }

        if ( opt->exists(key+".tolerance") ) {
            _epsTol = opt->get(key+".tolerance").as< double >();
        }
        else { _epsTol = 0.001; }    

}


bool XQMP::EvaluateFrame(Topology *top) {

    cout << endl << "... ... Register jobs, initialize MPS-mapper: " << flush;

    // CREATE XJOBS, INITIALIZE MAPPER
    _XJobs = XJOBS_FROM_TABLE(_job_file, top);
    _mps_mapper.GenerateMap(_xml_file, _emp_file, top, _XJobs);


    // RIGIDIFY TOPOLOGY (=> LOCAL FRAMES)
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else cout << endl << "... ... System is already rigidified." << flush;

    
    // CONVERT THREADS INTO SUBTHREADS IF BENEFICIAL
    if (_XJobs.size() < _nThreads) {
        _subthreads = (_nThreads - _XJobs.size()) / _XJobs.size() + 1;
        _nThreads   = _XJobs.size();

        cout << endl << "... ... "
             << "Converted threads into subthreads to increase efficiency: "
             << "NT = " << _nThreads << ", NST = " << _subthreads
             << flush;
    }

    
    // CREATE + EXECUTE THREADS (XJOB HANDLERS)
    vector<JobXQMP*> jobOps;
    _nextXJob = _XJobs.begin();

    for (int id = 0; id < _nThreads; id++) {
        JobXQMP *newOp = new JobXQMP(id, top, this);
        jobOps.push_back(newOp);
    }

    for (int id = 0; id < _nThreads; id++) {
        jobOps[id]->InitSlotData(top);
    }

    for (int id = 0; id < _nThreads; id++) {
        jobOps[id]->Start();
    }

    for (int id = 0; id < _nThreads; id++) {
        jobOps[id]->WaitDone();
    }

    for (int id = 0; id < _nThreads; id++) {
        delete jobOps[id];
    }

    jobOps.clear();

    
    // WRITE OUTPUT (PRIMARILY ENERGIE SPLITTINGS)
    FILE *out;
    out = fopen(this->_outFile.c_str(), "w");
    vector<XJob*> :: iterator jit;
    for (jit = _XJobs.begin(); jit < _XJobs.end(); ++jit) {
        (*jit)->WriteInfoLine(out);
    }
    fclose(out);
    
}


XJob *XQMP::RequestNextJob(int id, Topology *top) {

    _nextJobMutex.Lock();

    XJob *workOnThis;

    if (_nextXJob == _XJobs.end()) {
        workOnThis = NULL;
    }
    else {
        workOnThis = *_nextXJob;
        _nextXJob++;
        cout << endl 
             << "... ... Thread " << id << " evaluating job "
             << workOnThis->getId() << " " << workOnThis->getTag()
             << flush;
    }

    _nextJobMutex.Unlock();

    return workOnThis;
}


// ========================================================================== //
//                           JOBXQMP MEMBER FUNCTIONS                         //
// ========================================================================== //


void XQMP::JobXQMP::Run(void) {

    while (true) {
        _job = _master->RequestNextJob(_id, _top);

        if (_job == NULL) { break; }
        else { this->EvalJob(_top, _job); }
    }
}


void XQMP::JobXQMP::EvalJob(Topology *top, XJob *job) {
    
    // GENERATE POLAR TOPOLOGY
    double co1 = _master->_cutoff1;
    double co2 = _master->_cutoff2;    
    
    _master->_mps_mapper.Gen_QM_MM1_MM2(top, job, co1, co2);
    
    cout << endl << "... ... ... "
         << job->getPolarTop()->ShellInfoStr()
         << flush;
    
    if (tools::globals::verbose)
    job->getPolarTop()->PrintPDB(job->getTag()+"_QM0_MM1_MM2.pdb");

    // CALL MAGIC INDUCTOR         
    _inductor = XInductor(_master->_induce,     _master->_induce_intra_pair,
                          _master->_subthreads, _master->_wSOR_N,
                          _master->_wSOR_C,     _master->_epsTol,
                          _master->_maxIter,    _master->_maverick,
                          top,                  _id);
    
    _inductor.Evaluate(job);
    

    // SAVE INDUCTION STATE
    if (_master->_write_chk) {
        
        string format    = _master->_chk_format;
        string dotsuffix = (format == "gaussian") ? ".com" : ".xyz";
        string outstr    = job->getTag()+_master->_write_chk_suffix+dotsuffix;
        
        bool split       = _master->_chk_split_dpl;
        double space     = _master->_chk_dpl_spacing;
        
        job->getPolarTop()->PrintInduState(outstr, format, split, space);
        
    }

    // CLEAN POLAR TOPOLOGY    
    job->getPolarTop()->~PolarTop();
    
}


}}

#endif

