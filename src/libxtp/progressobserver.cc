#include <votca/xtp/progressobserver.h>
#include <votca/xtp/qmthread.h>
#include <boost/filesystem.hpp>
#include <boost/thread/thread.hpp>
#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

// Shuffling around #include directives between progobs.h and progobs.cc yields these errors:
// /people/thnfs/homes/poelking/VOTCA_SUSE_12/src/xtp/include/votca/xtp/logger.h:83:13: error: ‘string’ was not declared in this scope
// /usr/lib64/gcc/x86_64-suse-linux/4.7/../../../../x86_64-suse-linux/bin/ld: CMakeFiles/xtp_map.dir/xtp_map.cc.o: undefined reference to symbol _ZN5boost6system15system_categoryEv

using boost::format;

namespace votca { namespace xtp {
    
    
template<typename JobContainer, typename pJob, typename rJob>
pJob ProgObserver<JobContainer,pJob,rJob>::RequestNextJob(QMThread *thread) {
    
    _lockThread.Lock();    
    pJob jobToProc;
    
    LOG(logDEBUG,*(thread->getLogger())) 
        << "Requesting next job" << flush;

    // NEED NEW CHUNK?
    if (_nextjit == _jobsToProc.end() && _moreJobsAvailable) {
        SyncWithProgFile(thread);
        _nextjit = _jobsToProc.begin();
        if (_nextjit == _jobsToProc.end()) {
                 _moreJobsAvailable = false;
                LOG(logDEBUG,*(thread->getLogger()))
                        << "Sync did not yield any new jobs." << flush;
        }
    }
    
    // JOBS EATEN ALL UP?
    if (_nextjit == _jobsToProc.end()) {
        if (_maxJobs == _startJobsCount) {
            LOG(logDEBUG,*(thread->getLogger()))
                << "Next job: ID = - (reached maximum for this process)" 
                << flush;
        }
        else {
            LOG(logDEBUG,*(thread->getLogger())) 
                << "Next job: ID = - (none available)" << flush;
        }
        jobToProc = NULL;
    }
    // TAKE A BITE
    else {        
        jobToProc = *_nextjit;
        ++_nextjit;
        LOG(logDEBUG,*(thread->getLogger()))
            << "Next job: ID = " << jobToProc->getId() << flush;
    }
    
    if (!thread->isMaverick() && jobToProc != NULL) {
        int idx = jobToProc->getId();        
        int frac = (_jobs.size() >= 10) ? 10 : _jobs.size();
        int rounded = int(double(_jobs.size())/frac)*frac;
        int tenth = rounded / frac;        
        if (idx % tenth == 0) {
            double percent = double(idx-1) / rounded * 100 + 0.5;
            cout << (format("=> [%1$2.0f%%] ") % percent).str() << flush;
        }
    }
    
    _lockThread.Unlock();
    return jobToProc;
}


template<typename JobContainer, typename pJob, typename rJob>
void ProgObserver<JobContainer,pJob,rJob>::ReportJobDone(pJob job, rJob *res, QMThread *thread) {    
    _lockThread.Lock();
    LOG(logDEBUG,*(thread->getLogger()))
        << "Reporting job results" << flush;    
    // RESULTS, TIME, HOST
    job->SaveResults(res);    
    job->setTime(GenerateTime());
    job->setHost(GenerateHost(thread));
    // PRINT PROGRESS BAR
    _jobsReported += 1;
    if (!thread->isMaverick())
        cout << endl << *thread->getLogger() << flush;
    _lockThread.Unlock();
    return;
}


template<typename JobContainer, typename pJob, typename rJob>
string ProgObserver<JobContainer,pJob,rJob>::GenerateHost(QMThread *thread) {
    char host[128];
    //int h = 
    (void)gethostname(host, sizeof host);
    pid_t pid = getpid();
    //int tid = thread->getId(); // not used
    return (format("%1$s:%2$d") % host % pid).str();   
}


template<typename JobContainer, typename pJob, typename rJob>
string ProgObserver<JobContainer,pJob,rJob>::GenerateTime() {
    boost::posix_time::ptime now 
        = boost::posix_time::second_clock::local_time();
    return (format("%1$s") % now.time_of_day()).str();  
}


template<typename JobContainer, typename pJob, typename rJob>
void ProgObserver<JobContainer,pJob,rJob>::SyncWithProgFile(QMThread *thread) {
    
    // INTERPROCESS FILE LOCKING (THREAD LOCK IN ::RequestNextJob)
    this->LockProgFile(thread);
    
    string progFile = _progFile;
    string progBackFile = _progFile+"~";
    string tabFile = progFile;
    boost::algorithm::replace_last(tabFile, ".xml", ".tab");
    if (tabFile == progFile) tabFile += ".tab";
    string tabBackFile = tabFile+"~";
    
    // LOAD EXTERNAL JOBS FROM SHARED XML & UPDATE INTERNAL JOBS
    LOG(logDEBUG,*(thread->getLogger()))
        << "Update internal structures from job file" << flush;
    JobContainer jobs_ext = LOAD_JOBS<JobContainer,pJob,rJob>(progFile);    
    UPDATE_JOBS<JobContainer,pJob,rJob>(jobs_ext, _jobs, GenerateHost(thread));
    
    JobItVec it;
    for (it = jobs_ext.begin(); it != jobs_ext.end(); ++it) {
        pJob pj = *it;
        delete pj;
    }
    jobs_ext.clear();
    
    // GENERATE BACK-UP FOR SHARED XML
    LOG(logDEBUG,*(thread->getLogger()))
        << "Create job-file back-up" << flush;
    WRITE_JOBS<JobContainer,pJob,rJob>(_jobs, progBackFile, "xml");
    WRITE_JOBS<JobContainer,pJob,rJob>(_jobs, tabBackFile, "tab");
    
    
    // ASSIGN NEW JOBS IF AVAILABLE
    LOG(logDEBUG,*(thread->getLogger()))
        << "Assign jobs from stack" << flush;
    _jobsToProc.clear();
    
    unsigned int cacheSize = _cacheSize;
    while (_jobsToProc.size() < cacheSize) {
        if (_metajit == _jobs.end() || _startJobsCount == _maxJobs) break;
        
        bool startJob = false;
        
        // Start if job available or restart patterns matched
        if ( ((*_metajit)->isAvailable())
          || (_restartMode && _restart_stats.count((*_metajit)->getStatusStr()))
          || (_restartMode && _restart_hosts.count((*_metajit)->getHost())) )
            startJob = true;
        
        if (startJob) {
            (*_metajit)->Reset();
            (*_metajit)->setStatus("ASSIGNED");
            (*_metajit)->setHost(GenerateHost(thread));
            (*_metajit)->setTime(GenerateTime());
            _jobsToProc.push_back(*_metajit);
            _startJobsCount += 1;
        }

        ++_metajit;
    }
    
    // UPDATE PROGRESS STATUS FILE
    WRITE_JOBS<JobContainer,pJob,rJob>(_jobs, progFile, "xml");
    WRITE_JOBS<JobContainer,pJob,rJob>(_jobs, tabFile, "tab");

    // RELEASE PROGRESS STATUS FILE
    this->ReleaseProgFile(thread);
    return;
}


template<typename JobContainer, typename pJob, typename rJob>
void ProgObserver<JobContainer,pJob,rJob>::LockProgFile(QMThread *thread) {
    _flock = new boost::interprocess::file_lock(_lockFile.c_str());
    _flock->lock();
    LOG(logDEBUG,*(thread->getLogger()))
        << "Imposed lock on " << _lockFile << flush;
    LOG(logDEBUG,*(thread->getLogger()))
        << "Sleep ... " << _lockFile << flush;
    //boost::this_thread::sleep(boost::posix_time::milliseconds(0.0));
    LOG(logDEBUG,*(thread->getLogger()))
        << "Wake up ... " << _lockFile << flush;
}


template<typename JobContainer, typename pJob, typename rJob>
void ProgObserver<JobContainer,pJob,rJob>::ReleaseProgFile(QMThread *thread) {
    
    _flock->unlock();
    LOG(logDEBUG,*(thread->getLogger()))
        << "Releasing " << _lockFile << ". " << flush;
    delete _flock;
}


template<typename JobContainer, typename pJob, typename rJob>
void ProgObserver<JobContainer,pJob,rJob>
    ::InitCmdLineOpts(const boost::program_options::variables_map &optsMap) {
    
    _lockFile = optsMap["file"].as<string>();
    _cacheSize = optsMap["cache"].as<int>();
    _maxJobs = optsMap["maxjobs"].as<int>();
    string restartPattern = optsMap["restart"].as<string>();
    
    // restartPattern = e.g. host(pckr124:1234) stat(FAILED)    
    boost::algorithm::replace_all(restartPattern, " ", "");
    if (restartPattern == "") _restartMode = false;
    else _restartMode = true;
    
    vector<string> split;
    Tokenizer toker(restartPattern, "(,)");
    toker.ToVector(split);
    
    string category = "";
    for (unsigned int i = 0; i < split.size(); ++i) {
        
        if (split[i] == "host" || split[i] == "stat") category = split[i];
        
        else if (category == "host") _restart_hosts[split[i]] = true;
        else if (category == "stat") {
            if (split[i] == "ASSIGNED" || split[i] == "COMPLETE") 
                cout << "Restart if status == " << split[i] 
                    << "? Not necessarily a good idea." << endl;
            _restart_stats[split[i]] = true;
        }
        
        else throw runtime_error("Restart pattern ill-defined, format is"
                "[host([HOSTNAME:PID])] [stat([STATUS])]");
    }
    return;
}


template<typename JobContainer, typename pJob, typename rJob>
void ProgObserver<JobContainer,pJob,rJob>::InitFromProgFile(string progFile, 
    QMThread *thread) {
    
    _progFile = progFile;
    _jobsReported = 0;
    
    LOG(logINFO,*(thread->getLogger()))
        << "Job file = '" << _progFile << "', ";
    LOG(logINFO,*(thread->getLogger())) 
        << "lock file = '" << _lockFile << "', ";
    LOG(logINFO,*(thread->getLogger())) 
        << "cache size =  " << _cacheSize << flush;
    
    LOG(logINFO,*(thread->getLogger())) << "Initialize jobs from "
            << progFile << flush;    
    LOG(logINFO,*(thread->getLogger())) << "Lock & load " << flush;
    
    // LOCK, READ INTO XML
    this->LockProgFile(thread);  
    
    // ... Clear container
    JobItCnt it;
    for (it = _jobs.begin(); it != _jobs.end(); ++it) {
        pJob job = *it;
        delete job;
    }
    _jobs.clear();

    // ... Load new, set availability bool
    _jobs = LOAD_JOBS<JobContainer,pJob,rJob>(progFile);
    _metajit = _jobs.begin();
    WRITE_JOBS<JobContainer,pJob,rJob>(_jobs, progFile+"~", "xml");
    LOG(logINFO,*(thread->getLogger())) << "Registered " << _jobs.size()
         << " jobs." << flush;
	if (_jobs.size()>0) _moreJobsAvailable = true;
	else _moreJobsAvailable = false;
    
    
    // SUMMARIZE OBSERVER VARIABLES: RESTART PATTERN, CACHE, LOCK FILE
    if (_restartMode && _restart_hosts.size()) {        
        string infostr = "Restart if host == ";
        map<string,bool> ::iterator mit;
        for (mit = _restart_hosts.begin(); mit != _restart_hosts.end(); ++mit) {
            infostr += mit->first + " ";
        }
        LOG(logINFO,*(thread->getLogger())) << infostr << flush;  
    }
    if (_restartMode && _restart_stats.size()) {        
        string infostr = "Restart if stat == ";
        map<string,bool> ::iterator mit;
        for (mit = _restart_stats.begin(); mit != _restart_stats.end(); ++mit) {
            infostr += mit->first + " ";
        }
        LOG(logINFO,*(thread->getLogger())) << infostr << flush;  
    }   
    
    
    // RELEASE PROGRESS FILE
    this->ReleaseProgFile(thread);
    return;
}


// ========================================================================== //
//                         TEMPLATE SPECIALIZATIONS
// ========================================================================== //


template<typename JobContainer, typename pJob, typename rJob>
JobContainer LOAD_JOBS(const string &job_file) {    
    
    throw std::runtime_error("LOAD_JOBS not specialized for this type.");    
    JobContainer jobcnt;
    return jobcnt;
}

template<>
vector<Job*> LOAD_JOBS< vector<Job*>, Job*, Job::JobResult >(const string &job_file) {
    
    vector<Job*> jobs;    
    Property xml;
    load_property_from_xml(xml, job_file);
    
    list<Property*> jobProps = xml.Select("jobs.job");
    list<Property*> ::iterator it;
    for (it = jobProps.begin(); it != jobProps.end(); ++it) {
        
        Job *newJob = new Job(*it);
        jobs.push_back(newJob);       
    }
    
    return jobs;   
}

template<typename JobContainer, typename pJob, typename rJob>
void WRITE_JOBS(JobContainer &jobs, const string &job_file, string fileformat) {
    
    throw std::runtime_error("WRITE_JOBS not specialized for this type.");    
    return;
}

template<>
void WRITE_JOBS< vector<Job*>, Job*, Job::JobResult >(vector<Job*> &jobs, 
        const string &job_file, string fileformat) {
    
    vector<Job*> ::iterator it;
    
    ofstream ofs;    
    ofs.open(job_file.c_str(), ofstream::out);
    if (!ofs.is_open()) {
        throw runtime_error("Bad file handle: " + job_file);
    }    
    if (fileformat == "xml") ofs << "<jobs>" << endl;    
    for (it = jobs.begin(); it != jobs.end(); ++it) {
        if (fileformat == "tab" && !(*it)->isComplete()) continue;
        (*it)->ToStream(ofs, fileformat);
    }
    if (fileformat == "xml") ofs << "</jobs>" << endl;
    
    ofs.close();
    return;
}

template<typename JobContainer, typename pJob, typename rJob>
void UPDATE_JOBS(JobContainer &from, JobContainer &to, string thisHost) {
    
    throw std::runtime_error("UPDATE_JOBS not specialized for this type.");    
    return;
}

template<>
void UPDATE_JOBS< vector<Job*>, Job*, Job::JobResult >(vector<Job*> &from, 
        vector<Job*> &to, string thisHost) {
    
    vector<Job*> ::iterator it_int;
    vector<Job*> ::iterator it_ext;
    
    if (to.size() != from.size()) 
        throw runtime_error("Progress file out of sync (::size), abort.");
    
    for (it_int = to.begin(), it_ext = from.begin(); 
        it_int != to.end();
        ++it_int, ++it_ext) {
        
        Job* job_int = *it_int;
        Job* job_ext = *it_ext;
        
        if (job_int->getId() != job_ext->getId())
            throw runtime_error("Progress file out of sync (::id), abort.");
        
        if (job_ext->hasHost() && job_ext->getHost() != thisHost)
            job_int->UpdateFrom(job_ext);
    }
    
    return;
}




// REGISTER
template class ProgObserver< vector<Job*>, Job*, Job::JobResult >;
    
    
}}
