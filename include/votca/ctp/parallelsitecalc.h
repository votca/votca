#ifndef PARALLELSITECALC_H
#define PARALLELSITECALC_H


#include <votca/ctp/qmcalculator.h>
#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>


namespace votca { namespace ctp {

class ParallelSiteCalculator : public QMCalculator
{

public:

    ParallelSiteCalculator() : _nextSite(NULL) {};
   ~ParallelSiteCalculator() {};

    string       Identify() { return "Parallel Site Calculator"; }

    bool         EvaluateFrame(Topology *top);
    virtual void InitSlotData(Topology *top) { ; }
    virtual void PostProcess(Topology *top) { ; }
    virtual void EvalSite(Topology *top, Segment *seg, int slot) { ; }

    Segment     *RequestNextSite(int opId, Topology *top);
    void         LockCout() { _coutMutex.Lock(); }
    void         UnlockCout() { _coutMutex.Unlock(); }


    // ++++++++++++++++++++++++++++++++++++++ //
    // Site workers (i.e. individual threads) //
    // ++++++++++++++++++++++++++++++++++++++ //

    class SiteOperator : public Thread
    {
    public:

        SiteOperator(int id, Topology *top,
                     ParallelSiteCalculator *master)
                   : _id(id), _top(top), _seg(NULL),
                     _master(master)      {};

       ~SiteOperator() {};

        int  getId() { return _id; }
        void setId(int id) { _id = id; }

        void Run(void);
        

    protected:

        int                      _id;
        Topology                *_top;
        Segment                 *_seg;
        ParallelSiteCalculator  *_master;
    };


protected:

    vector<Segment*> ::iterator _nextSite;
    Mutex                       _nextSiteMutex;
    Mutex                       _coutMutex;
    

};

bool ParallelSiteCalculator::EvaluateFrame(Topology *top) {

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }
    cout << endl;    

    vector<SiteOperator*> siteOps;
    this->InitSlotData(top);

    _nextSite = top->Segments().begin();

    for (int id = 0; id < _nThreads; id++) {
        SiteOperator *newOp = new SiteOperator(id, top, this);
        siteOps.push_back(newOp);
    }

    for (int id = 0; id < _nThreads; id++) {
        siteOps[id]->Start();
    }

    for (int id = 0; id < _nThreads; id++) {
        siteOps[id]->WaitDone();
    }

    for (int id = 0; id < _nThreads; id++) {
        delete siteOps[id];
    }

    siteOps.clear();

    this->PostProcess(top);
    return 1;
}


// +++++++++++++++++ //
// Thread Management //
// +++++++++++++++++ //

Segment *ParallelSiteCalculator::RequestNextSite(int opId, Topology *top) {

    _nextSiteMutex.Lock();

    Segment *workOnThis;

    if (_nextSite == top->Segments().end()) {
        workOnThis = NULL;
    }
    else {
        workOnThis = *_nextSite;
        _nextSite++;
    }

    _nextSiteMutex.Unlock();

    return workOnThis;
}



// +++++++++++++++++++++++++++++ //
// SiteOperator Member Functions //
// +++++++++++++++++++++++++++++ //

void ParallelSiteCalculator::SiteOperator::Run(void) {

    while (true) {

        Segment *seg = _master->RequestNextSite(_id, _top);

        if (seg == NULL) { break; }
        else { this->_master->EvalSite(_top, seg, _id); }
    }
}


}}

#endif
