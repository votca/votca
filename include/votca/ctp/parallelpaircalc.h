#ifndef PARALLELPAIRCALC_H
#define PARALLELPAIRCALC_H


#include <votca/ctp/qmcalculator.h>
#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>


namespace votca { namespace ctp {

class ParallelPairCalculator : public QMCalculator
{

public:

    ParallelPairCalculator() : _nextPair(NULL) {};
   ~ParallelPairCalculator() {};

    string       Identify() { return "Parallel pair calculator"; }

    bool         EvaluateFrame(Topology *top);
    virtual void InitSlotData(Topology *top) { ; }
    virtual void PostProcess(Topology *top) { ; }
    virtual void EvalPair(Topology *top, QMPair *qmpair, int slot) { ; }

    QMPair     *RequestNextPair(int opId, Topology *top);
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
        QMPair                 *_pair;
        ParallelPairCalculator  *_master;
    };


protected:

    QMNBList::iterator   _nextPair;
    Mutex                 _nextPairMutex;
    Mutex                 _coutMutex;


};

bool ParallelPairCalculator::EvaluateFrame(Topology *top) {

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }
    cout << endl;        

    vector<PairOperator*> pairOps;
    this->InitSlotData(top);

    _nextPair = top->NBList().begin();

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

    this->PostProcess(top);
    return 1;
}


// +++++++++++++++++ //
// Thread Management //
// +++++++++++++++++ //

QMPair *ParallelPairCalculator::RequestNextPair(int opId, Topology *top) {

    _nextPairMutex.Lock();

    QMPair *workOnThis;

    if (_nextPair == top->NBList().end()) {
        workOnThis = NULL;
    }
    else {
        QMPair *workOnThat = *_nextPair;
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

        QMPair *qmpair = _master->RequestNextPair(_id, _top);

        if (qmpair == NULL) { break; }
        else { this->_master->EvalPair(_top, qmpair, _id); }
    }
}


}}





#endif
