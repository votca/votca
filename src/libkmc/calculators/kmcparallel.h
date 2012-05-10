#ifndef KMCPARALLEL_H
#define KMCPARALLEL_H


#include <votca/kmc/vssmgroup.h>
#include <vector>
#include <map>
#include <iostream>
#include <votca/tools/vec.h>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/globals.h>
#include "node.h"

#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>
#include <votca/tools/random2.h>


namespace votca { namespace kmc {
    
using namespace std;



class KMCParallel : public KMCCalculator
{
public:
    
    void        Initialize(const char *filename, Property *options );    
    bool        EvaluateFrame();
    bool        RequestNextInjection(int opId);
    
    
    // +++++++++++++++++++++ //
    // VSSM Group BOXED      //
    // +++++++++++++++++++++ //
    
    class NodeBoxed;
    
    template<typename event_t>
    class VSSMGroupBoxed
    {
    public:
        
        VSSMGroupBoxed() { _acc_rate.push_back(0); }
       ~VSSMGroupBoxed() {};

        void        AddEvent(event_t *event);
        double      Rate() { return _acc_rate.back(); };
        double      WaitingTime() { return _waiting_time; }
        void        UpdateWaitingTime();
        
        void        OnExecute();
        event_t    *getEvent(int idx) { return _events[idx]; }

    protected:

        event_t     *SelectEvent_LinearSearch();
        event_t     *SelectEvent_BinarySearch();

        std::vector<event_t*>       _events;
        std::vector<double>         _acc_rate;
        double                      _waiting_time;
        Random2                    *_random;
    };
    
    
    // +++++++++++++++++++++ //
    // Node BOXED            //
    // +++++++++++++++++++++ //
    
    struct LinkBoxed;
    
    class NodeBoxed : public VSSMGroupBoxed<LinkBoxed>
    {
    public:
        
        NodeBoxed(int id) : _id(id), _occ(0.0) {};
       ~NodeBoxed() {};
       
        void OnExecute() { 
            _occ += WaitingTime();
            VSSMGroupBoxed<LinkBoxed>::OnExecute();
        }
        
        int     _id;
        double  _occ;
    };   
    
    
    // +++++++++++++++++++++ //
    // Link BOXED            //
    // +++++++++++++++++++++ //
    
    class KMCSingleOp;
    
    struct LinkBoxed 
    {
        LinkBoxed(NodeBoxed *dest, double rate, vec dr, KMCSingleOp *master) 
                : _dest(dest), _rate(rate), _dr(dr), _master(master) {};
       ~LinkBoxed() {};
       
       double     Rate() { return _rate; }
       void       OnExecute() { _master->Step(_dr, _dest); }
       
       double           _rate;
       NodeBoxed       *_dest;
       vec              _dr;
       KMCSingleOp     *_master;
    };
    
    
    // +++++++++++++++++++++ //
    // KMC OPERATOR          //
    // +++++++++++++++++++++ //
    
    class KMCSingleOp : public Thread
    {
    public:
        
        KMCSingleOp(int id, KMCParallel *master) : _id(id), _master(master) {};
       ~KMCSingleOp() {};
       
        int     getId() { return _id; }
        void    setId(int id) { _id = id; }
        
        void    InitSlotData();
        void    Run(void);
        
        void    EvalKMC();
        void    LoadGraph();
        void    RunKMC(void);
        void    WriteOcc(void);
        void    Reset();
        
        void    Step(vec &dr, NodeBoxed *dest) { _pos += dr; _current = dest; }
        
    
    private:
        
        int                     _id;
        KMCParallel            *_master; 
        
        map<int , node_t *>     _nodes_lookup;
        vector<node_t *>        _nodes;
        map<int , node_t *>     _injection_lookup;
        vector<node_t *>        _injection;  
        
        vec                     _pos;
        NodeBoxed              *_current;
    };
    
    
private:    

    string      _stateFile;   
   
    // KMC run variables
    string      _injName;
    string      _channel;
    double      _runtime;
    double      _outtime;
    int         _NRuns;
    int         _nextRun;
    int         _seed;
    Random2     _seeder;
    
    // KMC log variables
    string      _outFile;
    bool        _kmc2File;
    map<int, vector<vec> > _log_inj_vel;
    
    // Thread management
    Mutex       _injMutex;
    Mutex       _logMutex;
    bool        _maverick;
    
};


// ++++++++++++++++++++ //
// KMC PARALLEL MEMBERS //
// ++++++++++++++++++++ //


void KMCParallel::Initialize(const char* filename, Property *options) {
    
    cout << "Initialize "
         << flush;
    
    cout << endl
         << "... KMCParallel "
         << flush;   
    
    cout << endl 
         << "... ... Initializing KMCParallel with " << _nThreads << " threads."
         << flush;
    
    // Init. master parameters
    string key  = "options.kmcparallel";
    _stateFile  = filename;
    _maverick   = (_nThreads == 1) ? true : false;

    _runtime    = options->get(key+".runtime").as<double>();
    _outtime    = options->get(key+".outtime").as<double>();
    
    _seed       = options->get(key+".seed").as<int>();
    _injName    = options->get(key+".injection").as<string>();
    _channel    = options->get(key+".channel").as<string>();
    _NRuns      = options->get(key+".runs").as<int>();
    _nextRun    = 1;
    
    // Init. master RNG
    srand(_seed);
    _seeder.init(rand(), rand(), rand(), rand());
    
}


bool KMCParallel::EvaluateFrame() {
    
    cout << "Evaluate Frame " << flush;
    cout << endl << "... KMCParallel" << flush;
    
    vector<KMCSingleOp*> kmcOps;
    
    for (int id = 0; id < _nThreads; ++id) {
        KMCSingleOp *newOp = new KMCSingleOp(id, this);
        kmcOps.push_back(newOp);
    }
    
    for (int id = 0; id < _nThreads; ++id) {
        kmcOps[id]->InitSlotData();
    }
    
    for (int id = 0; id < _nThreads; ++id) {
        kmcOps[id]->Start();
    }
    
    for (int id = 0; id < _nThreads; ++id) {
        kmcOps[id]->WaitDone();
    }
    
    for (int id = 0; id < _nThreads; ++id) {
        delete kmcOps[id];
    }
    
    kmcOps.clear();
    
    cout << endl;
    
}

bool KMCParallel::RequestNextInjection(int opId) {
    
    _injMutex.Lock();
    
    bool doNext = false;
    
    if (_nextRun > _NRuns) {
        doNext = false;
    }
    else {
        
        cout << "\r"
             << "... ... OP " << opId << " "
             << "starting injection run " 
             << _nextRun << "/" << _NRuns
             << ". " << flush;
        
        ++_nextRun;
        doNext = true;
    }
    
    _injMutex.Unlock();
    
    return doNext;
}


    
// ++++++++++++++++++++ //
// KMC OPERATOR MEMBERS //
// ++++++++++++++++++++ //    


void KMCParallel::KMCSingleOp::InitSlotData() {
    ;
}


void KMCParallel::KMCSingleOp::LoadGraph() {
    ;
}


void KMCParallel::KMCSingleOp::Run(void) {
    
    while (true) {
        
        bool doNext = _master->RequestNextInjection(this->_id);
        
        if (!doNext) { break; }
        else { this->EvalKMC(); }        
    }
}


void KMCParallel::KMCSingleOp::EvalKMC() {    
    ;
}


void KMCParallel::KMCSingleOp::WriteOcc() {
    ;
}


void KMCParallel::KMCSingleOp::Reset() {
    ;
}


// ++++++++++++++++++++ //
// VSSM GROUP MEMBERS   //
// ++++++++++++++++++++ //

template<typename event_t>
void KMCParallel::VSSMGroupBoxed<event_t>::AddEvent(event_t *event) {
    
    _events.push_back(event);
    _acc_rate.push_back(_acc_rate.back() + event->Rate());
    
    UpdateWaitingTime();
}


template<typename event_t>
inline void KMCParallel::VSSMGroupBoxed<event_t>::UpdateWaitingTime() {
    
    _waiting_time = -log( 1.0 - _random->rand_uniform() ) / Rate();
}
  

template<typename event_t>
inline void KMCParallel::VSSMGroupBoxed<event_t>::OnExecute() {
    
    SelectEvent_BinarySearch()->OnExecute();
    UpdateWaitingTime();
}


template<typename event_t>
event_t *KMCParallel::VSSMGroupBoxed<event_t>::SelectEvent_BinarySearch() {
        
    double max = Rate();
    double u = 1.-_random->rand_uniform();
    u=u*Rate();	

    int imin = 0;
    int imax = _acc_rate.size();

    while(imax - imin > 1) {
            int imid = (int) ((imin+imax)*0.5);
            if(u <= _acc_rate[imid]) imax=imid;
            else                     imin=imid;
    }

    return _events[imin];
}

}}

#endif
