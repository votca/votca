/*
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#ifndef PARALLELSITECALC_H
#define PARALLELSITECALC_H


#include <votca/xtp/qmcalculator.h>
#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>


namespace votca { namespace xtp {

class ParallelSiteCalculator : public QMCalculator
{

        
public:

    class SiteOperator;
    ParallelSiteCalculator() : _nextSite(NULL) {};
   ~ParallelSiteCalculator() {};

    string       Identify() { return "Parallel Site Calculator"; }

    bool         EvaluateFrame(Topology *top);
    virtual void InitSlotData(Topology *top) { ; }
    virtual void PostProcess(Topology *top) { ; }
    virtual void EvalSite(Topology *top, Segment *seg, int slot, SiteOperator* opThread ) { ; }

    Segment     *RequestNextSite(int opId, Topology *top);
    void         LockCout() { _coutMutex.Lock(); }
    void         UnlockCout() { _coutMutex.Unlock(); }


    // ++++++++++++++++++++++++++++++++++++++ //
    // Site workers (i.e. individual threads) //
    // ++++++++++++++++++++++++++++++++++++++ //

    class SiteOperator : public QMThread
    {
    public:

        SiteOperator(int id, Topology *top,
                     ParallelSiteCalculator *master)
                   : _top(top), _seg(NULL),
                     _master(master)      { _id = id; };

       ~SiteOperator() {};

        void Run(void);
        

    protected:

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
        //cout << endl << "... ... " << "Evaluating site " << workOnThis->getId() << endl;
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
        else { this->_master->EvalSite(_top, seg, _id, this); }
    }
}


}}

#endif
