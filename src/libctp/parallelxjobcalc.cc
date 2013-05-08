#include <votca/ctp/parallelxjobcalc.h>
#include <boost/algorithm/string.hpp>


namespace votca { namespace ctp {

    
template<typename JobContainer, typename pJob> 
bool ParallelXJobCalc<JobContainer,pJob>::EvaluateFrame(Topology *top) {    

    // CREATE XJOBS & PROGRESS OBSERVER (_XJOBFILE INIT. IN CHILD)
    _XJobs = XJOBS_FROM_TABLE<JobContainer,pJob>(_xjobfile, top); 
    cout << endl << "... ... Registered " << _XJobs.size() << " jobs " << flush;
    
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
    
    string progFile = _xjobfile+"_status";
    string lockFile = _xjobfile;
    _progObs = new ProgObserver< JobContainer, pJob >
            (&_XJobs, _nThreads, progFile, lockFile);


    // PRE-PROCESS (OVERWRITTEN IN CHILD OBJECT)
    this->PreProcess(top);
    
    // CREATE + EXECUTE THREADS (XJOB HANDLERS)
    vector<JobOperator*> jobOps;

    for (int id = 0; id < _nThreads; id++) {
        JobOperator *newOp = new JobOperator(id, top, this);
        jobOps.push_back(newOp);
    }
    
    for (int id = 0; id < _nThreads; ++id) {
        CustomizeLogger(jobOps[id]);
    }

    for (int id = 0; id < _nThreads; id++) {
        jobOps[id]->InitData(top);
    }

    for (int id = 0; id < _nThreads; id++) {
        jobOps[id]->Start();
    }

    for (int id = 0; id < _nThreads; id++) {
        jobOps[id]->WaitDone();
    }
    
    if (!_maverick)
    for (int id = 0; id < _nThreads; id++) {
        cout << endl << *(jobOps[id]->getLogger()) << flush;
    }

    for (int id = 0; id < _nThreads; id++) {
        delete jobOps[id];
    }    

    jobOps.clear();

    
    // POST-PROCESS (OVERWRITTEN IN CHILD OBJECT)
    this->PostProcess(top);
    
}


template<typename JobContainer, typename pJob>
void ParallelXJobCalc<JobContainer,pJob>::JobOperator::Run(void) {

    while (true) {
        _job = _master->_progObs->RequestNextJob(this);

        if (_job == NULL) { break; }
        else { 
            this->_master->EvalJob(_top, _job, this); 
            this->_master->_progObs->ReportJobDone(_job, this);
        }
    }
}


// REGISTER PARALLEL CALCULATORS
template class ParallelXJobCalc< vector<XJob*>, XJob* >;
template class ParallelXJobCalc< vector<Segment*>, Segment* >;

}}