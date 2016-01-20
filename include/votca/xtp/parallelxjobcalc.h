#ifndef __PARALLELXJOBCALC__H
#define __PARALLELXJOBCALC__H


#include <votca/xtp/jobcalculator.h>
#include <votca/xtp/qmthread.h>
#include <votca/tools/mutex.h>
#include <votca/xtp/job.h>
#include <votca/xtp/progressobserver.h>


// PATHWAYS TO A NEW THREADED CALCULATOR
// ... 1 Define 'JobContainer' (needs to define iterator), 'pJob' ( = *iterator)
// ... 2 Derive new calculator as ': public ParallelXJobCalc<JobContainer,pJob>'
// ... 3 Specialize XJOBS_FROM_TABLE< JobContainer, pJob> in xjob.cc
// ... 4 Register new calculator (see end of parallelxjobcalc.cc)

// REQUIRED MEMBERS FOR pJob
// pJob::JobResult (struct)


namespace votca { namespace xtp {

template<typename JobContainer, typename pJob, typename rJob> 
class ParallelXJobCalc : public JobCalculator
{

public:

    class JobOperator;
    
    ParallelXJobCalc() : _jobfile("__NOFILE__") {};
   ~ParallelXJobCalc() { ; };

    string       Identify() { return "Parallel XJob Calculator"; }

    bool         EvaluateFrame(Topology *top);
    virtual void LoadJobs() { ; }
    virtual void CustomizeLogger(QMThread* thread);
    virtual void PreProcess(Topology *top) { ; } 
    virtual rJob EvalJob(Topology *top, const pJob job, QMThread *thread) = 0;
    virtual void PostProcess(Topology *top) { ; }
    
    void         LockCout() { _coutMutex.Lock(); }
    void         UnlockCout() { _coutMutex.Unlock(); }
    void         LockLog() { _logMutex.Lock(); }
    void         UnlockLog() { _logMutex.Unlock(); }

    
    // ======================================== //
    // XJOB OPERATOR (THREAD)                   //
    // ======================================== //
    

    class JobOperator : public QMThread
    {
    public:

        JobOperator(int id,   Topology *top, ParallelXJobCalc<JobContainer,pJob,rJob> *master)
                      : _top(top),          _master(master) { _id = id; };
       ~JobOperator() {};

        void        InitData(Topology *top) { ; }
        void        Run(void);
        

    public:

        Topology         *_top;
        ParallelXJobCalc<JobContainer,pJob,rJob> *_master;
        pJob              _job;

    };

    
    
    


protected:

    JobContainer             _XJobs;
    Mutex                    _coutMutex;
    Mutex                    _logMutex;
    string                   _jobfile;
    int                      _subthreads;
    
    // ProgObserver< JobContainer, pJob > *_progObs;


};

}}

#endif