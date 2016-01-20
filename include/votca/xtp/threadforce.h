#ifndef VOTCA_XTP_THREADFORCE_H
#define VOTCA_XTP_THREADFORCE_H

#include <votca/tools/thread.h>

namespace votca {
namespace xtp {

template
<
    // REQUIRED Clone()
    typename T
>
class PrototypeCreator
{
public:
    PrototypeCreator(T *prototype=0) : _prototype(prototype) {}
    T *Create() { return (_prototype) ? _prototype->Clone() : 0; }
    T *getPrototype() { return _prototype; }
    void setPrototype(T *prototype) { _prototype = prototype; }
private:
    T *_prototype;
};


template
<
    // REQUIRED Start(), WaitDone()
    // OPTIONAL AddAtomicInput(<>), AddSharedInput(<>),
    //          setMode(<>), setVerbose(bool)
    typename ThreadType,
    // REQUIRED Create()
    template<typename> class CreatorType
>
class ThreadForce : 
    public vector<ThreadType*>, 
    public CreatorType<ThreadType>
{
public:
    typedef typename ThreadForce::iterator it_t;
    
    ThreadForce() {}
   ~ThreadForce() { Disband(); }
    
    void Initialize(int n_threads) {
        // Set-up threads via CreatorType policy
        for (int t = 0; t < n_threads; ++t) {
            ThreadType *new_thread = this->Create();
            //new_thread->setId(this->size()+1);
            this->push_back(new_thread);
        }
        return;
    }
    
    template<typename mode_t>
    void AssignMode(mode_t mode) {
        // Set mode
        it_t tfit;
        for (tfit = this->begin(); tfit < this->end(); ++tfit)
           (*tfit)->setMode(mode);
        return;
    }
    
    template<class in_t>
    void AddAtomicInput(vector<in_t> &inputs) {
        // Apportion the vector elements equally onto all threads
        typename vector<in_t>::iterator iit;
        int in_idx = 0;
        int th_idx = 0;
        for (iit = inputs.begin(); iit != inputs.end(); ++iit, ++in_idx) {
            th_idx = in_idx % this->size();
            (*this)[th_idx]->AddAtomicInput(*iit);
        }
        // The thread with the last input can be verbose
        for (it_t tfit = this->begin(); tfit < this->end(); ++tfit)
            (*tfit)->setVerbose(false);
        (*this)[th_idx]->setVerbose(true);
        return;
    }
    
    template<class in_t>
    void AddSharedInput(in_t &shared_input) {
        for (it_t tfit = this->begin(); tfit < this->end(); ++tfit)
            (*tfit)->AddSharedInput(shared_input);
        return;
    }
    
    void StartAndWait() {
        // Start threads
        for (it_t tfit = this->begin(); tfit < this->end(); ++tfit)
            (*tfit)->Start();
        // Wait for threads
        for (it_t tfit = this->begin(); tfit < this->end(); ++tfit)
            (*tfit)->WaitDone();
        return;
    }
    
    template<class mode_t>
    void Reset(mode_t mode) {
        for (it_t tfit = this->begin(); tfit < this->end(); ++tfit)
            (*tfit)->Reset(mode);
        return;
    }
    
    void Disband() {
        // Delete & clear
        it_t tfit;
        for (tfit = this->begin(); tfit < this->end(); ++tfit)
            delete *tfit;
        this->clear();
        return;
    }
    
private:
    
};
    

template
<
    // REQUIRED void Run(void)
    typename ThreadBase,
    // REQUIRED (none)
    typename ThreadDerived
>
class MultiModeThread : public ThreadBase
{
public:
    
    MultiModeThread() {};
   ~MultiModeThread() {};
   
    typedef void (MultiModeThread::*StartFunct)();
    typedef void (MultiModeThread::*ResetFunct)();
    typedef double (MultiModeThread::*WloadFunct)();
   
    const string &getMode() { return _current_mode; }
    void setMode(const string &mode) { _current_mode = mode; }
    void setVerbose(bool verbose) { _verbose = verbose; }
   
    void RegisterStart(string mode, void (ThreadDerived::*member)()) {
        _mode_startfunct[mode] = static_cast<StartFunct>(member);
    }
    
    void RegisterReset(string mode, void (ThreadDerived::*member)()) {
        _mode_resetfunct[mode] = static_cast<ResetFunct>(member);
    }
    
    void RegisterWload(string mode, double (ThreadDerived::*member)()) {
        _mode_wloadfunct[mode] = static_cast<WloadFunct>(member);
    }    
    
    void Run(void) {
        StartFunct start = _mode_startfunct[_current_mode]; 
        ((*this).*start)();
    }
   
    void Reset(string mode) {
        ResetFunct reset = _mode_resetfunct[mode];
        ((*this).*reset)();
    }
    
    double Workload(string mode) {
        WloadFunct wload = _mode_wloadfunct[mode];
        return ((*this).*wload)();
    }
    
protected:
    
    bool _verbose;
    string _current_mode;
    map<string,StartFunct> _mode_startfunct;
    map<string,ResetFunct> _mode_resetfunct;
    map<string,WloadFunct> _mode_wloadfunct;
        
};
    
    
}}




#endif