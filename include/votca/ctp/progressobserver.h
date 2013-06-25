#ifndef VOTCA_CTP_PROGRESSOBSERVER
#define VOTCA_CTP_PROGRESSOBSERVER


#include <vector>
#include <iostream>
#include <votca/tools/mutex.h>
#include <votca/tools/property.h>
#include <votca/ctp/job.h>
#include <boost/interprocess/sync/file_lock.hpp>

using namespace std;

namespace votca { namespace ctp {
    
class QMThread;

// TYPENAME EXAMPLE USAGE
//     ProgObserver< vector<Job*>, Job*, Job::JobResult >
// REQUIRED METHODS FOR TYPENAMES
//     pJob ->getId() ->SaveResults(rJob)
//     JobContainer .size() .begin() .end()
    
template<typename JobContainer, typename pJob, typename rJob>
class ProgObserver 
{
    
public:
    
    typedef typename JobContainer::iterator JobItCnt;
    typedef typename vector<pJob>::iterator JobItVec;
    
    ProgObserver(int nThreads, string stateFile)
        : _nThreads(nThreads), _lockFile(stateFile) { ; }
    
    ProgObserver()
        : _nThreads(-1), _progFile("__NOFILE__"), 
          _lockFile("__NOFILE__"),    _nextjit(NULL), _metajit(NULL) { ; }
    
   ~ProgObserver() { ; }
    
    void InitFromProgFile(string progFile, QMThread *master);   
    pJob RequestNextJob(QMThread *thread);
    void ReportJobDone(pJob job, rJob *res, QMThread *thread);
    
    void SyncWithProgFile(QMThread *thread);
    void LockProgFile(QMThread *thread);
    void ReleaseProgFile(QMThread *thread);
    
    string GenerateStamp(QMThread *thread);
    string GenerateTime();
   
   
private:    
    
    JobItCnt _metajit;
    JobContainer _jobs;
    
    JobItVec _nextjit;
    vector<pJob> _jobsToProc;
    vector<pJob> _jobsToSync;
    
    int _nThreads;
    string _progFile;
    string _lockFile;
    Mutex _lockThread;
    boost::interprocess::file_lock *_flock;
    
};




template<typename JobContainer, typename pJob, typename rJob>
JobContainer LOAD_JOBS(const string &xml_file);

template<typename JobContainer, typename pJob, typename rJob>
void WRITE_JOBS(JobContainer &jobs, const string &job_file, string fileformat);

template<typename JobContainer, typename pJob, typename rJob>
void UPDATE_JOBS(JobContainer &from, JobContainer &to);
    
    
    
    
}}


#endif
