#ifndef KMCPARALLEL_H
#define KMCPARALLEL_H


#include <vector>
#include <map>
#include <iostream>
#include <votca/tools/vec.h>
#include <votca/tools/database.h>
#include <votca/tools/statement.h>
#include <votca/tools/tokenizer.h>
#include <votca/tools/globals.h>

#include <votca/tools/thread.h>
#include <votca/tools/mutex.h>
#include <votca/tools/random2.h>


namespace votca { namespace xtp {
    
using namespace std;



class KMCParallel : public KMCCalculator
{
public:
    
    using       KMCCalculator::Initialize;
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

        void        SetRNG(Random2 *rng) { _random = rng; }
       
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
        
        void Reset() { _occ = 0.0; }
        
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
       
       NodeBoxed       *_dest;
       double           _rate;
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
        Random2                 _random;
        
        map<int,NodeBoxed* >    _nodes_lookup;
        vector< NodeBoxed* >    _nodes;
        map<int,NodeBoxed* >    _injection_lookup;
        vector< NodeBoxed* >    _injection;  
        
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
         << "... ... Initialize KMCParallel with " << _nThreads << " threads. "
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
    
    _outFile    = options->get(key+".output").as<string>();
    
    // Init. master RNG
    srand(_seed);
    cout << endl << "... ... Master: " << flush;
    _seeder.init(rand(), rand(), rand(), rand());
    
}


bool KMCParallel::EvaluateFrame() {
    
    cout << "Evaluate Frame " << endl;
    cout << "... KMCParallel" << endl;
    
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
    
    
    FILE *out;
    out = fopen(_outFile.c_str(),"w");
    
    map< int, vector<vec> > ::iterator mit;
    vector<vec> ::iterator vit;
    vec AVGVEL = vec(0,0,0);
    int AVG_N  = 0;
    
    for (mit = _log_inj_vel.begin();
         mit != _log_inj_vel.end();
         ++mit) {
        for (vit = (*mit).second.begin(); 
             vit < (*mit).second.end();
             ++vit) {
            fprintf(out, "%5d %4.7e %4.7e %4.7e \n",
              (*mit).first, (*vit).getX(), (*vit).getY(), (*vit).getZ());
            AVGVEL += (*vit);
            AVG_N  += 1;
        }
    }   
    
    AVGVEL /= AVG_N;
    fprintf(out, "AVG VELOCITY %4.7e %4.7e %4.7e \n", 
            AVGVEL.getX(), AVGVEL.getY(), AVGVEL.getZ());    
    fclose(out);
    
    cout << "... ... AVG VELOCITY = " << AVGVEL << " m/s." << endl;
    return true;
}


bool KMCParallel::RequestNextInjection(int opId) {
    
    _injMutex.Lock();
    
    bool doNext = false;
    
    if (_nextRun > _NRuns) {
        doNext = false;
    }
    else {
        
        cout << ""
             << "... ... OP " << opId << " "
             << "starting injection run " 
             << _nextRun << "/" << _NRuns
             << ". " << endl;
        
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
    
    // Initialise random-number generator
    cout << "... ... OP " << this->_id << ": " << flush;
    this->_random.init(rand(), rand(), rand(), rand());
    
    Database db;
    db.Open(_master->_stateFile);
    
    if (_master->_maverick) {
        cout << flush
             << "... ... Loading graph from " << _master->_stateFile
             << endl;
    }
    
    // Load segments <=> nodes
    Statement *stmt = db.Prepare("SELECT id, name FROM segments;");
    
    while (stmt->Step() != SQLITE_DONE) {
        
        int     id      = stmt->Column<int>(0);
        string  name    = stmt->Column<string>(1);
        
        NodeBoxed *newnode = new NodeBoxed(id);
        newnode->SetRNG(&this->_random);
        
        _nodes.push_back(newnode);
        _nodes_lookup[id] = _nodes.back();
        
        if (wildcmp(_master->_injName.c_str(), name.c_str())) {
            _injection.push_back(newnode);
            _injection_lookup[id] = _injection.back();
        }
    }
    
    delete stmt;
    
    // Load pairs <=> links
    int linkCount = 0;
    if (_master->_channel == "hole") {
        stmt = db.Prepare("SELECT seg1,    seg2, "
                                 "rate12h, rate21h, "
                                 "drX, drY, drZ "
                          "FROM   pairs;");
    }
    else if (_master->_channel == "electron") {
        stmt = db.Prepare("SELECT seg1,    seg2, "
                                 "rate12e, rate21e, "
                                 "drX, drY, drZ "
                          "FROM   pairs;");
    }
    else {
        throw std::runtime_error(" Invalid channel option '" +
                                   _master->_channel + "'. ");
    }
    
    while (stmt->Step() != SQLITE_DONE) {
        
        NodeBoxed *node1 = _nodes_lookup[stmt->Column<int>(0)];
        NodeBoxed *node2 = _nodes_lookup[stmt->Column<int>(1)];
        
        double rate12 = stmt->Column<double>(2);
        double rate21 = stmt->Column<double>(3);
        
        vec dr = vec(stmt->Column<double>(4), 
                     stmt->Column<double>(5),
                     stmt->Column<double>(6));
        
        node1->AddEvent(new LinkBoxed(node2, rate12, dr, this));
        node2->AddEvent(new LinkBoxed(node1, rate21, -dr, this));
        ++(++linkCount);
    }
    
    delete stmt;
    stmt = NULL;
    
    if (_master->_maverick) {
        cout << flush
             << "... ... Created graph with "
             << _nodes.size() << " nodes, "
             << linkCount << " links. "
             << endl;
    }    
}


void KMCParallel::KMCSingleOp::LoadGraph() {
    ; // ::InitSlotData()
}


void KMCParallel::KMCSingleOp::Run(void) {
    
    while (true) {
        
        bool doNext = _master->RequestNextInjection(this->_id);
        
        if (!doNext) { break; }
        else { this->EvalKMC(); }        
    }
}


void KMCParallel::KMCSingleOp::EvalKMC() {

    // Pick injection site
    int inj = _random.rand_uniform_int(_injection.size());
    _current = _injection[inj];
    _current->UpdateWaitingTime();
    _pos = vec(0,0,0);
    
    if (_master->_maverick) {
        cout << "... ... ... Starting run at node " 
             << _current->_id << " (" << inj << "). "
             << endl;
    }
    
    // Init KMC
    int    steps = 0;
    double t_run = 0.0;
    double t_out = _master->_outtime;
    double t_max = _master->_runtime;
        
    // Prepare trajectory log file
    FILE *out;
    string logTrajFile = "log_traj_" 
                       + boost::lexical_cast<string>(_id) 
                       + ".dat";
    out = fopen(logTrajFile.c_str(), "a");
    fprintf(out, "Injection at node %5d \n", inj);
    fprintf(out, "%5d %4.7e %4.7f %4.7f %4.7f \n", 
            _current->_id, t_run, _pos.getX(), _pos.getY(), _pos.getZ());
    
    // Run KMC
    while (t_run < t_max) {
        
        t_run += _current->WaitingTime();
        _current->OnExecute();
        
        if (t_run > t_out) {
            t_out = t_run + _master->_outtime;
            fprintf(out, "%5d %4.7e %4.7f %4.7f %4.7f \n", 
            _current->_id, t_run, _pos.getX(), _pos.getY(), _pos.getZ());            
        }
        
        ++steps;
    }
    fprintf(out, "%5d %4.7e %4.7f %4.7f %4.7f \n", 
    _current->_id, t_run, _pos.getX(), _pos.getY(), _pos.getZ()); 
    fprintf(out, "Finished KMC in %9d steps. \n", steps);
    fprintf(out, "AVG velocity %4.7e %4.7e %4.7e \n\n",
                   (_pos/t_run*1e-9).getX(), 
                   (_pos/t_run*1e-9).getY(), 
                   (_pos/t_run*1e-9).getZ());
    fclose(out);
    
    _master->_logMutex.Lock();
    _master->_log_inj_vel[inj].push_back(_pos/t_run*1e-9);
    _master->_logMutex.Unlock();
    
    // Post-process
    this->WriteOcc();
    this->Reset();
}


void KMCParallel::KMCSingleOp::WriteOcc() {
    ;
}


void KMCParallel::KMCSingleOp::Reset() {
    vector< NodeBoxed* > ::iterator nit;
    for (nit = _nodes.begin(); nit < _nodes.end(); ++nit) {
        (*nit)->Reset();
    }
    
    _pos = vec(0.,0.,0.);
    _current = NULL;
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
        
    //double max = Rate();
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
