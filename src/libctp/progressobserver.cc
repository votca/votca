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

using boost::format;

namespace votca { namespace ctp {
    

template<typename JobContainer, typename pJob>
pJob ProgObserver<JobContainer,pJob>::RequestNextJob(QMThread *thread) {
    
    _lockRequest.Lock();    
    pJob jobToProc;
    
    LOG(logDEBUG,*(thread->getLogger())) << "Requesting next job" << flush;
    
    // NEED NEW CHUNK?
    if (_nextjit == _jobsToProc.end()) {
        SyncWithProgFile(thread);
        _nextjit = _jobsToProc.begin();
    }
    
    // JOBS EATEN ALL UP?
    if (_nextjit == _jobsToProc.end()) {
        LOG(logDEBUG,*(thread->getLogger())) << "Request: No more." << flush;
        jobToProc = NULL;
    }
    // TAKE A BITE
    else {        
        jobToProc = *_nextjit;
        ++_nextjit;
        LOG(logDEBUG,*(thread->getLogger())) << "Request, so more: " << jobToProc->getId() << flush;
    }
    
    _lockRequest.Unlock();
    return jobToProc;
}


template<typename JobContainer, typename pJob>
void ProgObserver<JobContainer,pJob>::SyncWithProgFile(QMThread *thread) {
    
    // INTERPROCESS FILE LOCKING (THREAD LOCKING IN ::RequestNextJob)
    this->LockProgFile(thread);
   
    string progFile = _progFile;
    ofstream ofs;
    ifstream ifs;
    
    // WRITE PROGRESS STATUS FILE IF NECESSARY
    if (!boost::filesystem::exists(progFile)) {
        ofs.open(progFile.c_str(), ofstream::out);
        if (!ofs.is_open()) {
            LOG(logERROR,*(thread->getLogger())) << "Could not open file "
                << progFile << flush;
            throw runtime_error("Bad file handle.");
        }
        LOG(logDEBUG,*(thread->getLogger())) << "Created new file "
            << progFile << flush;
        JobItCnt jit;
        for (jit = _jobs->begin(); jit != _jobs->end(); ++jit) {
            ofs << WriteProgLine(*jit, thread, "QUEUEING") << endl;
        }
        ofs.close();
    }
    
    // READ PROGRESS STATUS FILE INTO MAPS
    map<int,string> assigned;
    map<int,string> complete;
    map<int,string> queueing;
    map<int,string> ::iterator atJobId;
    
    ifs.open(progFile.c_str(), ifstream::in);    
    if (!ifs.is_open()) {
        LOG(logERROR,*(thread->getLogger())) << "Could not open file "
            << progFile << flush;
        throw runtime_error("Bad file handle.");
    }
    while (ifs.good()) {
        string line;
        std::getline(ifs, line);
        vector< string > split;
        Tokenizer toker(line, " \t");
        toker.ToVector(split);
        if ( !split.size()      ||
              split[0] == "#"   ||
              split[0].substr(0,1) == "#" ) { continue; }

        int id = boost::lexical_cast<int>(split[0]);
        string status = split[1];
        
        if (status == "COMPLETE")
            complete[id] = line;
        else if (status == "ASSIGNED")
            assigned[id] = line;
        else if (status == "QUEUEING")
            queueing[id] = line;
        else
            assert(false); // unrecognized status        
    }    
    ifs.close();    
    
    // REPORT FINISHED JOBS, CLEAR SYNC CONTAINER
    LOG(logDEBUG,*(thread->getLogger())) << "Sync finished jobs " << flush;
    JobItVec vjit;
    for (vjit = _jobsToSync.begin(); vjit != _jobsToSync.end(); ++vjit) {
        pJob job = *vjit;
        atJobId = assigned.find(job->getId());
        if (atJobId == assigned.end())
            throw runtime_error("Job recently finished, but had not been assigned.");
        complete[atJobId->first] = WriteProgLine(job, thread, "COMPLETE");
        assigned.erase(atJobId);
    }
    _jobsToSync.clear();
    
    // BACK-UP PROGRESS FILE
    string backupfile = _progFile+"_backup";
    //if (boost::filesystem::exists(backupfile))
    //    boost::filesystem::remove(backupfile);
    //boost::filesystem::copy_file(_progFile, backupfile);
    LOG(logDEBUG,*(thread->getLogger())) 
        << "Back-up progress file for the contingency of a restart ..." << flush;
    ofs.open(backupfile.c_str(), ofstream::out);
    ofs << "# Verify that no jobs are missing (e.g. cat ... | wc) "
        << "before restarting from this file." << endl;
    JobItCnt jit;
    for (atJobId = complete.begin(); atJobId != complete.end(); ++atJobId) {
        pJob job = (*_jobs)[atJobId->first-1];
        ofs << atJobId->second << endl;        
    }
    for (atJobId = assigned.begin(); atJobId != assigned.end(); ++atJobId) {
        pJob job = (*_jobs)[atJobId->first-1];
        string modStr = atJobId->second;
        boost::replace_first(modStr, "ASSIGNED", "QUEUEING");
        ofs << modStr << endl;        
    }
    for (atJobId = queueing.begin(); atJobId != queueing.end(); ++atJobId) {
        pJob job = (*_jobs)[atJobId->first-1];
        ofs << atJobId->second << endl;        
    }
    ofs.close();    
    
    // ASSIGN NEW JOBS IF AVAILABLE
    LOG(logDEBUG,*(thread->getLogger())) << "Assign jobs from stack" << flush;
    _jobsToProc.clear();
    while (_jobsToProc.size() < _nThreads) {
        if (_metajit == _jobs->end()) break;        
        
        atJobId = queueing.find((*_metajit)->getId());
        
        if (atJobId == queueing.end()) {            
            ; // Job already processed elsewhere
        }
        else {
            _jobsToProc.push_back(*_metajit);
            assigned[atJobId->first] = WriteProgLine(*_metajit, thread, "ASSIGNED");
            queueing.erase(atJobId);
        }
        
        ++_metajit;
    }
    
    // SANITY CHECKS
    int jobCountMaps = queueing.size() + assigned.size() + complete.size();
    assert(jobCountMaps == _jobs->size());        
    
    // UPDATE PROGRESS STATUS FILE
    LOG(logDEBUG,*(thread->getLogger())) << "Update progress file ..." << flush;
    ofs.open(progFile.c_str(), ofstream::out);
    for (atJobId = complete.begin(); atJobId != complete.end(); ++atJobId) {
        pJob job = (*_jobs)[atJobId->first-1];
        ofs << atJobId->second << endl;        
    }
    for (atJobId = assigned.begin(); atJobId != assigned.end(); ++atJobId) {
        pJob job = (*_jobs)[atJobId->first-1];
        ofs << atJobId->second << endl;        
    }
    for (atJobId = queueing.begin(); atJobId != queueing.end(); ++atJobId) {
        pJob job = (*_jobs)[atJobId->first-1];
        ofs << atJobId->second << endl;        
    }
    ofs.close();

    // RELEASE PROGRESS STATUS FILE
    this->ReleaseProgFile(thread);
    
//    // REPORT FINISHED JOBS, CLEAR SYNC CONTAINER
//    string progFileFinished = _progFile + "_finished";    
//    
//    // Create file if not yet there
//    if (!boost::filesystem::exists(progFileFinished)) {
//        ofs.open(progFileFinished.c_str(), ofstream::out);
//        LOG(logERROR,*(thread->getLogger())) << "Created new file "
//            << progFileFinished << flush;
//    }
//    else {
//        ofs.open(progFileFinished.c_str(), ofstream::out | ofstream::app);
//    }
//    
//    if (!ofs.is_open()) {
//        LOG(logERROR,*(thread->getLogger())) << "Could not open file "
//            << progFileFinished << flush;
//        throw runtime_error("Bad file handle.");
//    }
//    
//    JobItVec vjit;
//    for (vjit = _jobsToSync.begin(); vjit != _jobsToSync.end(); ++vjit) {
//        ofs << WriteProgLine(*vjit, "finished");
//    }
//    
//    ofs.close();    
//    _jobsToSync.clear();
//    
//    
//    // READ PREVIOUSLY ASSIGNED JOBS INTO MAP
//    string progFileAssigned = _progFile + "_assigned";
//    
//    // Create file if not yet there
//    if (!boost::filesystem::exists(progFileAssigned)) {
//        ofs.open(progFileAssigned.c_str(), ofstream::out);
//        ofs.close();
//        LOG(logERROR,*(thread->getLogger())) << "Created new file "
//            << progFileAssigned << flush;
//    }
//    
//    map<int, string> assignedId_line;
//
//    ifs.open(progFileAssigned.c_str(), ifstream::in);
//
//    if (!ifs.is_open()) {
//        LOG(logERROR,*(thread->getLogger())) << "Could not open file "
//            << progFileAssigned << flush;
//        throw runtime_error("Bad file handle.");
//    }    
//
//    while ( ifs.good() ) {
//        string line;
//        std::getline(ifs, line);
//
//        vector< string > split;
//        Tokenizer toker(line, " \t");
//        toker.ToVector(split);
//        if ( !split.size()      ||
//              split[0] == "#"   ||
//              split[0].substr(0,1) == "#" ) { continue; }
//
//        int id = boost::lexical_cast<int>(split[0]);
//        assignedId_line[id] = line;
//    }
//    
//    ifs.close();
//    
//    // GRAB NEW JOBS TO PROCESS
//    _jobsToProc.clear();
//    
//    while(_jobsToProc.size() < _nThreads) {
//        
//        if (_metajit == _jobs->end()) break;
//        
//        
//        if ( assignedId_line.count((*_metajit)->getId()) > 0 ) {
//            ; // Already being processed somewhere
//        }
//        else {
//            _jobsToProc.push_back(*_metajit);
//        }
//        
//        ++_metajit;        
//    }
//    
//    // APPEND NEWLY ASSIGNED JOBS TO LIST OF ASSIGNED JOBS
//    ofs.open(progFileAssigned.c_str(), ofstream::out | ofstream::app);   
//    
//    if (!ofs.is_open()) {
//        LOG(logERROR,*(thread->getLogger())) << "Could not open file "
//            << progFileAssigned << flush;
//        throw runtime_error("Bad file handle.");
//    }
//    
//    for (vjit = _jobsToProc.begin(); vjit != _jobsToProc.end(); ++vjit) {
//        ofs << WriteProgLine(*vjit, "assigned");
//    }
//    
//    ofs.close();    
    
}


template<typename JobContainer, typename pJob>
void ProgObserver<JobContainer,pJob>::LockProgFile(QMThread *thread) {
    
    _flock = new boost::interprocess::file_lock(_lockFile.c_str());
    _flock->lock();
    LOG(logDEBUG,*(thread->getLogger()))
        << "Imposed lock on " << _lockFile << flush;
}


template<typename JobContainer, typename pJob>
void ProgObserver<JobContainer,pJob>::ReleaseProgFile(QMThread *thread) {
    
    _flock->unlock();
    LOG(logDEBUG,*(thread->getLogger()))
        << "Releasing " << _lockFile << ". " << flush;
    delete _flock;
}


template<typename JobContainer, typename pJob>
string ProgObserver<JobContainer,pJob>::WriteProgLine(pJob job, 
    QMThread *thread, string status) {
    
    pid_t pid = getpid();
    int   tid = thread->getId();
    
    string jobstr0 = "";
    string jobstr1 = "";
    string jobstr2 = "";
    
    // HOST NAME
    char host[128];
    int h = gethostname(host, sizeof host);
    
    // CURRENT TIME
    boost::posix_time::ptime now 
        = boost::posix_time::second_clock::local_time();
    
    // Job ID + Status
    jobstr0 += (format("%1$5d %2$8s ")
        % job->getId() % status).str();
    
    if (status == "QUEUEING") {
        jobstr1 += (format("HOST=........ TIME=........ PID=..... TID=... ")).str();
    }
    
    else if (status == "ASSIGNED") {
        jobstr1 += (format("HOST=%2$-8s TIME=%3$s PID=%1$05d TID=... ")
            % pid % host % now.time_of_day()).str();
        // += job-specific 'assigned' line
    }
    
    else if (status == "COMPLETE") {
        jobstr1 += (format("HOST=%3$-8s TIME=%4$s PID=%1$05d TID=%2$03d ")
        % pid % tid % host % now.time_of_day()).str();
        // += jobstr2 -> job-specific output line
    }
    
    else ;
    
    return jobstr0 + jobstr1 + jobstr2;
}


template<typename JobContainer, typename pJob>
void ProgObserver<JobContainer,pJob>::ReportJobDone(pJob job, QMThread *thread) {
    
    LOG(logDEBUG,*(thread->getLogger())) << "Reporting job done." << flush;
    _jobsToSync.push_back(job);
}


// ========================================================================== //
//                         TEMPLATE SPECIALIZATIONS
// ========================================================================== //


template<>
string ProgObserver< vector<XJob*>, XJob* >::WriteProgLine(XJob *job, 
    QMThread *thread, string status) {
    
    pid_t pid = getpid();
    int   tid = thread->getId();
    
    string jobstr0 = "";
    string jobstr1 = "";
    string jobstr2 = "";
    
    // HOST NAME
    char host[128];
    int h = gethostname(host, sizeof host);
    
    // CURRENT TIME
    boost::posix_time::ptime now 
        = boost::posix_time::second_clock::local_time();
    
    // Job ID + Status
    jobstr0 += (format("%1$5d %2$8s ")
        % job->getId() % status).str();
    
    if (status == "QUEUEING") {
        jobstr1 += (format("HOST=........ TIME=........ PID=..... TID=... ")).str();
    }
    
    else if (status == "ASSIGNED") {
        jobstr1 += (format("HOST=%2$-8s TIME=%3$s PID=%1$05d TID=... ")
            % pid % host % now.time_of_day()).str();
        jobstr2 += (format("%1$5d %2$-10s ")
            % job->getUserId() % job->getTag()).str();
    }
    
    else if (status == "COMPLETE") {
        jobstr1 += (format("HOST=%3$-8s TIME=%4$s PID=%1$05d TID=%2$03d ")
        % pid % tid % host % now.time_of_day()).str();
        jobstr2 += job->getInfoLine();
    }
    
    else ;
    
    return jobstr0 + jobstr1 + jobstr2;
}


// REGISTER
template class ProgObserver< vector<XJob*>, XJob* >;
template class ProgObserver< vector<Segment*>, Segment* >;
    
    
    
    
    
    
}}
