/*
 *            Copyright 2009-2018 The VOTCA Development Team
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
/// For an earlier history see ctp repo commit 77795ea591b29e664153f9404c8655ba28dc14e9


#ifndef VOTCA_CTP_PROGRESSOBSERVER
#define VOTCA_CTP_PROGRESSOBSERVER


#include <vector>
#include <iostream>
#include <votca/tools/mutex.h>
#include <votca/tools/property.h>
#include <votca/ctp/job.h>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/program_options.hpp>

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
    
    ProgObserver()
        : _lockFile("__NOFILE__"), _progFile("__NOFILE__"), _cacheSize(-1),
          _nextjit(NULL), _metajit(NULL), _startJobsCount(0) { ; }
    
   ~ProgObserver() { ; }
    
    void InitCmdLineOpts(const boost::program_options::variables_map &optsMap);
    void InitFromProgFile(string progFile, QMThread *master);   
    pJob RequestNextJob(QMThread *thread);
    void ReportJobDone(pJob job, rJob *res, QMThread *thread);
    
    void SyncWithProgFile(QMThread *thread);
    void LockProgFile(QMThread *thread);
    void ReleaseProgFile(QMThread *thread);
    
    string GenerateHost(QMThread *thread);
    string GenerateTime();
   
   
private:    
    
    string _lockFile;
    string _progFile;
    int _cacheSize;
    JobContainer _jobs;
    
    vector<pJob> _jobsToProc;
    vector<pJob> _jobsToSync;
    
    JobItVec _nextjit;
    JobItCnt _metajit;
    Mutex _lockThread;
    boost::interprocess::file_lock *_flock;
    
    map<string,bool> _restart_hosts;
    map<string,bool> _restart_stats;
    bool _restartMode;
    int _jobsReported;

    bool _moreJobsAvailable;
    int _startJobsCount;
    int _maxJobs;
    
};




template<typename JobContainer, typename pJob, typename rJob>
JobContainer LOAD_JOBS(const string &xml_file);

template<typename JobContainer, typename pJob, typename rJob>
void WRITE_JOBS(JobContainer &jobs, const string &job_file, string fileformat);

template<typename JobContainer, typename pJob, typename rJob>
void UPDATE_JOBS(JobContainer &from, JobContainer &to, string thisHost);
    
    
    
    
}}


#endif
