#ifndef PARALLELPAIRCALC_H
#define PARALLELPAIRCALC_H


#include <votca/ctp/qmcalculator2.h>
#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>


namespace votca { namespace ctp {

class ParallelPairCalculator : public QMCalculator2
{

public:

    ParallelPairCalculator() : _nThreads(1), _nextPair(NULL) {};
   ~ParallelPairCalculator() {};

    string       Identify() { return "Parallel pair calculator"; }

    bool         EvaluateFrame(Topology *top);
    virtual void PrepareFrame(Topology *top) { ; }
    virtual void FinishFrame(Topology *top) { ; }
    virtual void EvalPair(Topology *top, QMPair2 *qmpair, int slot) { ; }

    QMPair2     *RequestNextPair(int opId, Topology *top);
    void         LockCout() { _coutMutex.Lock(); }
    void         UnlockCout() { _coutMutex.Unlock(); }


    // ++++++++++++++++++++++++++++++++++++++ //
    // Pair workers (i.e. individual threads) //
    // ++++++++++++++++++++++++++++++++++++++ //

    class PairOperator : public Thread
    {
    public:

        PairOperator(int id, Topology *top,
                     ParallelPairCalculator *master)
                   : _id(id), _top(top), _pair(NULL),
                     _master(master)      {};

       ~PairOperator() {};

        int  getId() { return _id; }
        void setId(int id) { _id = id; }

        void Run(void);
        

    protected:

        int                      _id;
        Topology                *_top;
        QMPair2                 *_pair;
        ParallelPairCalculator  *_master;
    };


protected:

    int                   _nThreads;
    QMNBList2::iterator   _nextPair;
    Mutex                 _nextPairMutex;
    Mutex                 _coutMutex;


};

bool ParallelPairCalculator::EvaluateFrame(Topology *top) {

    _nextPair = top->NBList().begin();

    this->PrepareFrame(top);
    cout << endl;

    vector<PairOperator*> pairOps;

    for (int id = 0; id < _nThreads; id++) {
        PairOperator *newOp = new PairOperator(id, top, this);
        pairOps.push_back(newOp);
    }

    for (int id = 0; id < _nThreads; id++) {
        pairOps[id]->Start();
    }

    for (int id = 0; id < _nThreads; id++) {
        pairOps[id]->WaitDone();
    }

    for (int id = 0; id < _nThreads; id++) {
        delete pairOps[id];
    }

    pairOps.clear();

    this->FinishFrame(top);
    return 1;
}


// +++++++++++++++++ //
// Thread Management //
// +++++++++++++++++ //

QMPair2 *ParallelPairCalculator::RequestNextPair(int opId, Topology *top) {

    _nextPairMutex.Lock();

    QMPair2 *workOnThis;

    if (_nextPair == top->NBList().end()) {
        workOnThis = NULL;
    }
    else {
        QMPair2 *workOnThat = *_nextPair;
        _nextPair++;
        workOnThis = workOnThat;
    }

    _nextPairMutex.Unlock();

    return workOnThis;
}

// +++++++++++++++++++++++++++++ //
// PairOperator Member Functions //
// +++++++++++++++++++++++++++++ //

void ParallelPairCalculator::PairOperator::Run(void) {

    while (true) {

        QMPair2 *qmpair = _master->RequestNextPair(_id, _top);

        if (qmpair == NULL) { break; }
        else { this->_master->EvalPair(_top, qmpair, _id); }
    }
}


}}





#endif
