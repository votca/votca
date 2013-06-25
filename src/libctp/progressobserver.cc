#include <votca/ctp/progressobserver.h>
#include <votca/ctp/qmthread.h>
#include <boost/filesystem.hpp>
#include <boost/thread/thread.hpp>
#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

// Shuffling around #include directives between progobs.h and progobs.cc yields these errors:
// /people/thnfs/homes/poelking/VOTCA_SUSE_12/src/ctp/include/votca/ctp/logger.h:83:13: error: ‘string’ was not declared in this scope
// /usr/lib64/gcc/x86_64-suse-linux/4.7/../../../../x86_64-suse-linux/bin/ld: CMakeFiles/ctp_map.dir/ctp_map.cc.o: undefined reference to symbol _ZN5boost6system15system_categoryEv

using boost::format;

namespace votca { namespace ctp {
    
    
template<typename JobContainer, typename pJob, typename rJob>
pJob ProgObserver<JobContainer,pJob,rJob>::RequestNextJob(QMThread *thread) {
    
    _lockThread.Lock();    
    pJob jobToProc;
    
    LOG(logDEBUG,*(thread->getLogger())) 
        << "Requesting next job" << flush;
    
    // NEED NEW CHUNK?
    if (_nextjit == _jobsToProc.end()) {
        SyncWithProgFile(thread);
        _nextjit = _jobsToProc.begin();
    }
    
    // JOBS EATEN ALL UP?
    if (_nextjit == _jobsToProc.end()) {
        LOG(logDEBUG,*(thread->getLogger())) 
            << "Next job: ID = - (none available)" << flush;
        jobToProc = NULL;
    }
    // TAKE A BITE
    else {        
        jobToProc = *_nextjit;
        ++_nextjit;
        LOG(logDEBUG,*(thread->getLogger()))
            << "Next job: ID = " << jobToProc->getId() << flush;
    }
    
    _lockThread.Unlock();
    return jobToProc;
}


template<typename JobContainer, typename pJob, typename rJob>
void ProgObserver<JobContainer,pJob,rJob>::ReportJobDone(pJob job, rJob *res, QMThread *thread) {    
    _lockThread.Lock();    
    LOG(logDEBUG,*(thread->getLogger()))
        << "Reporting job results" << flush;    
    // RESULTS, TIME, STAMP
    job->SaveResults(res);    
    job->setTime(GenerateTime());
    job->setStamp(GenerateStamp(thread));
    _lockThread.Unlock();
    return;
}


template<typename JobContainer, typename pJob, typename rJob>
string ProgObserver<JobContainer,pJob,rJob>::GenerateStamp(QMThread *thread) {
    char host[128];
    int h = gethostname(host, sizeof host);
    pid_t pid = getpid();
    int tid = thread->getId();
    return (format("%1$s:%2$d:%3$d") % host % pid % tid).str();   
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
    UPDATE_JOBS<JobContainer,pJob,rJob>(jobs_ext, _jobs);
    
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
    while (_jobsToProc.size() < _nThreads) {
        if (_metajit == _jobs.end()) break;        
        
        if ((*_metajit)->isAvailable()) {            
            (*_metajit)->setStatus("ASSIGNED");
            (*_metajit)->setStamp(GenerateStamp(thread));
            (*_metajit)->setTime(GenerateTime());
            _jobsToProc.push_back(*_metajit);
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
}


template<typename JobContainer, typename pJob, typename rJob>
void ProgObserver<JobContainer,pJob,rJob>::ReleaseProgFile(QMThread *thread) {
    
    _flock->unlock();
    LOG(logDEBUG,*(thread->getLogger()))
        << "Releasing " << _lockFile << ". " << flush;
    delete _flock;
}


template<typename JobContainer, typename pJob, typename rJob>
void ProgObserver<JobContainer,pJob,rJob>::InitFromProgFile(string progFile, 
    QMThread *thread) {
    
    _progFile = progFile;
    
    LOG(logINFO,*(thread->getLogger())) << "Initialize jobs from "
            << progFile << flush;    
    LOG(logINFO,*(thread->getLogger())) << "Lock & load " << flush;
    
    // LOCK, READ INTO XML
    this->LockProgFile(thread);  
    
    // TODO Delete jobs here before loading new
    _jobs = LOAD_JOBS<JobContainer,pJob,rJob>(progFile);
    _metajit = _jobs.begin();
    WRITE_JOBS<JobContainer,pJob,rJob>(_jobs, progFile+"~", "xml");
    LOG(logINFO,*(thread->getLogger())) << "Registered " << _jobs.size()
         << " jobs." << flush;
    
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
    for (it = jobs.begin(); it != jobs.end(); ++it)  
        (*it)->ToStream(ofs, fileformat);
    if (fileformat == "xml") ofs << "</jobs>" << endl;
    
    ofs.close();
    return;
}

template<typename JobContainer, typename pJob, typename rJob>
void UPDATE_JOBS(JobContainer &from, JobContainer &to) {
    
    throw std::runtime_error("UPDATE_JOBS not specialized for this type.");    
    return;
}

template<>
void UPDATE_JOBS< vector<Job*>, Job*, Job::JobResult >(vector<Job*> &from, vector<Job*> &to) {
    
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
        
        job_int->UpdateFrom(job_ext);
    }
    
    return;
}




// REGISTER
template class ProgObserver< vector<Job*>, Job*, Job::JobResult >;
    
    
}}
