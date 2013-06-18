#ifndef VOTCA_CTP_PROGRESSOBSERVER
#define VOTCA_CTP_PROGRESSOBSERVER


#include <vector>
#include <iostream>
#include <votca/tools/mutex.h>
#include <votca/tools/property.h>
#include <votca/ctp/job.h>
#include <votca/ctp/xjob.h>
#include <boost/interprocess/sync/file_lock.hpp>

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
    
    ProgObserver(int nThreads, string stateFile)
        : _nThreads(nThreads), _lockFile(stateFile) { ; }
    
    ProgObserver()
        : _jobs(NULL), _nThreads(-1), _progFile("nofile"), _lockFile("nofile"),
          _nextjit(NULL), _metajit(NULL) { ; }
    
    ProgObserver(JobContainer *jobs, int nThreads, string sharedProgFile, string lockFile)
        : _jobs(jobs), _nThreads(nThreads), _progFile(sharedProgFile), _lockFile(lockFile)
          { _metajit = _jobs->begin(); _nextjit = _jobsToProc.begin(); }
    
   ~ProgObserver() { ; }
   
   string getLockFile() { return _lockFile; }
    void InitFromProgFile(string progFile, QMThread *master);
   
    pJob RequestNextJob(QMThread *thread);
    void ReportJobDone(pJob job, QMThread *thread);
    
    void SyncWithProgFile(QMThread *thread);
    void LockProgFile(QMThread *thread);
    string WriteProgLine(pJob job, QMThread *thread, string status);
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
    string _lockFile;
    Mutex _lockRequest;
    boost::interprocess::file_lock *_flock;
    
}; 
    
    
    
    
}}


#endif
