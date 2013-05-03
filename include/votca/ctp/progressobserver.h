#ifndef VOTCA_CTP_PROGRESSOBSERVER
#define VOTCA_CTP_PROGRESSOBSERVER


#include <vector>
#include <iostream>
#include <votca/tools/mutex.h>
#include <votca/ctp/xjob.h>

using namespace std;

namespace votca { namespace ctp {

class QMThread;
    
// E.G. ProgObserver< vector<Seg*>, Seg* >
//      ProgObserver< vector<XJob*>, XJob* >
//      ProgObserver< QMNBList, QMPair* >    
    
template<typename JobContainer, typename pJob>
class ProgObserver 
{
    
public:
    
    typedef typename JobContainer::iterator JobItCnt;
    typedef typename vector<pJob>::iterator JobItVec;
    
    ProgObserver()
        : _jobs(NULL), _nThreads(-1), _progFile("nofile"),
          _nextjit(NULL), _metajit(NULL) { ; }
    
    ProgObserver(JobContainer *jobs, int nThreads, string sharedProgFile)
        : _jobs(jobs), _nThreads(nThreads), _progFile(sharedProgFile)
          { _metajit = _jobs->begin(); _nextjit = _jobsToProc.begin(); }
    
   ~ProgObserver() { ; }
   
   
    pJob RequestNextJob(QMThread *thread);
    void ReportJobDone(pJob job, QMThread *thread);
    
    void SyncWithProgFile(QMThread *thread);
    void LockProgFile(QMThread *thread);
    string WriteProgLine(pJob job, string progStage);
    void ReleaseProgFile(QMThread *thread);
    void ReportJobOutcome(pJob job, QMThread *thread) {;}
   
   
private:    
    
    JobItCnt _metajit;
    JobContainer *_jobs;
    
    JobItVec _nextjit;
    vector<pJob> _jobsToProc;
    vector<pJob> _jobsToSync;
    
    int _nThreads;
    string _progFile;
    Mutex _lockRequest;
    
}; 
    
    
    
    
}}


#endif
