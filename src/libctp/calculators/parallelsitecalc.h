#ifndef PARALLELSITECALC_H
#define PARALLELSITECALC_H


#include <votca/ctp/qmcalculator2.h>
#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>


namespace votca { namespace ctp {

class ParallelSiteCalculator : public QMCalculator2
{

public:

    ParallelSiteCalculator() : _nThreads(1), _nextSite(1) { };
   ~ParallelSiteCalculator() { };

    string Identify() { return "Parallel Site Calculator"; }

    bool     EvaluateFrame(Topology *top);
    Segment *RequestNextSite(int opId, Topology *top);
    void     LockCout() { _coutMutex.Lock(); }
    void     UnlockCout() { _coutMutex.Unlock(); }


    // +++++++++++++++++++++++++++++++++ //
    // Threaded workers (site operators) //
    // +++++++++++++++++++++++++++++++++ //

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
        virtual void EvalSite(Topology *top, Segment *seg);

    protected:

        int                      _id;
        Topology                *_top;
        Segment                 *_seg;
        ParallelSiteCalculator  *_master;
    };


protected:

    int       _nThreads;
    int       _nextSite;
    Mutex     _nextSiteMutex;
    Mutex     _coutMutex;


};


bool ParallelSiteCalculator::EvaluateFrame(Topology *top) {

    cout << endl;

    vector<SiteOperator*> siteOps;

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

    return 1;
}


// +++++++++++++++++ //
// Thread Management //
// +++++++++++++++++ //

Segment *ParallelSiteCalculator::RequestNextSite(int opId, Topology *top) {

    _nextSiteMutex.Lock();

    Segment *workOnThis;

    if (_nextSite < 1) {        
        workOnThis = NULL;
    }
    else {
        workOnThis = top->getSegment(_nextSite);
        _nextSite++;
        if (_nextSite > top->Segments().size()) {
            _nextSite = -1;
        }
    }

    _nextSiteMutex.Unlock();

    return workOnThis;
}



// +++++++++++++++++++++++++++++ //
// SiteOperator Member Functions //
// +++++++++++++++++++++++++++++ //

void ParallelSiteCalculator::SiteOperator::Run(void) {

    while (true) {

        Segment *seg = _master->RequestNextSite(this->getId(), _top);

        if (seg == NULL) { break; }
        else { this->EvalSite(_top, seg); }
    }
}

/*
virtual void ParallelSiteCalculator::SiteOperator::EvalSite(Topology *top, Segment *seg)
{

    this->_master->LockCout();
    cout << "\r... ... Evaluating site " << seg->getId() << ". " << flush;
    this->_master->UnlockCout();

    int ij;
    for (int i = 0; i < 2000; i++) {
        for (int j = 0; j < 2000; j++) {
            ij = i+j;
        }
    }
}
*/



}}

#endif
