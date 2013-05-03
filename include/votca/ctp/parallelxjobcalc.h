#ifndef __PARALLELXJOBCALC__H
#define __PARALLELXJOBCALC__H


#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/qmthread.h>
#include <votca/tools/mutex.h>
#include <votca/ctp/xjob.h>
#include <votca/ctp/progressobserver.h>


namespace votca { namespace ctp {

class ParallelXJobCalc : public QMCalculator
{

public:

    class XJobOperator;
    
    ParallelXJobCalc() : _nextXJob(NULL), _xjobfile("__NOFILE__") {};
   ~ParallelXJobCalc() {};

    string       Identify() { return "Parallel XJob Calculator"; }

    bool         EvaluateFrame(Topology *top);
    virtual void CustomizeLogger(XJobOperator* thread) { ; }
    virtual void PreProcess(Topology *top) { ; } 
    virtual void EvalJob(Topology *top, XJob *qmpair, XJobOperator* opThread) { ; }
    virtual void PostProcess(Topology *top) { ; }
    
    void         LockCout() { _coutMutex.Lock(); }
    void         UnlockCout() { _coutMutex.Unlock(); }
    void         LockLog() { _logMutex.Lock(); }
    void         UnlockLog() { _logMutex.Unlock(); }

    
    // ======================================== //
    // XJOB OPERATOR (THREAD)                   //
    // ======================================== //
    

    class XJobOperator : public QMThread
    {
    public:

        XJobOperator(int id,   Topology *top, ParallelXJobCalc *master)
                      : _id(id),    _top(top),          _master(master) {};
       ~XJobOperator() {};

        int         getId() { return _id; }
        void        setId(int id) { _id = id; }
        void        InitData(Topology *top) { ; }
        void        Run(void);
        

    public:

        int               _id;
        Topology         *_top;
        ParallelXJobCalc *_master;
        XJob             *_job;

    };

    
    
    


protected:

    vector<XJob*>            _XJobs;
    vector<XJob*> ::iterator _nextXJob;
    Mutex                    _nextJobMutex;
    Mutex                    _coutMutex;
    Mutex                    _logMutex;
    bool                     _maverick;
    string                   _xjobfile;
    int                      _subthreads;
    
    ProgObserver< vector<XJob*>, XJob* > _progObs;


};

}}

#endif