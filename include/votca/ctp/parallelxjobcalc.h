#ifndef __PARALLELXJOBCALC__H
#define __PARALLELXJOBCALC__H


#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/qmthread.h>
#include <votca/tools/mutex.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/progressobserver.h>


namespace votca { namespace ctp {

template<typename JobContainer, typename pJob> 
class ParallelXJobCalc : public QMCalculator
{

public:

    class JobOperator;
    
    ParallelXJobCalc() : _xjobfile("__NOFILE__") {};
   ~ParallelXJobCalc() { delete _progObs; };

    string       Identify() { return "Parallel XJob Calculator"; }

    bool         EvaluateFrame(Topology *top);
    virtual void CustomizeLogger(QMThread* thread) { ; }
    virtual void PreProcess(Topology *top) { ; } 
    virtual void EvalJob(Topology *top, XJob *qmpair, QMThread *opThread) { ; }
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

        JobOperator(int id,   Topology *top, ParallelXJobCalc<JobContainer,pJob> *master)
                      : _top(top),          _master(master) { _id = id; };
       ~JobOperator() {};

        void        InitData(Topology *top) { ; }
        void        Run(void);
        

    public:

        Topology         *_top;
        ParallelXJobCalc<JobContainer,pJob> *_master;
        pJob              _job;

    };

    
    
    


protected:

    vector<XJob*>            _XJobs;
    Mutex                    _coutMutex;
    Mutex                    _logMutex;
    bool                     _maverick;
    string                   _xjobfile;
    int                      _subthreads;
    
    ProgObserver< JobContainer, pJob > *_progObs;


};

}}

#endif