#include <votca/ctp/parallelxjobcalc.h>


namespace votca { namespace ctp {


bool ParallelXJobCalc::EvaluateFrame(Topology *top) {

    // CREATE XJOBS FROM FILE (HAS TO BE INITIALIZED IN CHILD OBJECT)
    _XJobs = XJOBS_FROM_TABLE(_xjobfile, top);    
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

    // PRE-PROCESS (OVERWRITTEN IN CHILD OBJECT)
    this->PreProcess(top);
    
    // CREATE + EXECUTE THREADS (XJOB HANDLERS)
    vector<XJobOperator*> jobOps;
    _nextXJob = _XJobs.begin();

    for (int id = 0; id < _nThreads; id++) {
        XJobOperator *newOp = new XJobOperator(id, top, this);
        jobOps.push_back(newOp);
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

    for (int id = 0; id < _nThreads; id++) {
        delete jobOps[id];
    }

    jobOps.clear();

    
    // POST-PROCESS (OVERWRITTEN IN CHILD OBJECT)
    this->PostProcess(top);
    
}


XJob *ParallelXJobCalc::RequestNextJob(int id, Topology *top) {

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


void ParallelXJobCalc::XJobOperator::Run(void) {

    while (true) {
        _job = _master->RequestNextJob(_id, _top);

        if (_job == NULL) { break; }
        else { this->_master->EvalJob(_top, _job, this); }
    }
}

}}