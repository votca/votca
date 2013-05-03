#include <votca/ctp/progressobserver.h>
#include <votca/ctp/qmthread.h>
#include <boost/filesystem.hpp>
#include <boost/thread/thread.hpp>
#include <boost/format.hpp>
#include <fstream>

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
    
    this->LockProgFile(thread);
   
    
    // REPORT FINISHED JOBS, CLEAR SYNC CONTAINER
    string progFileFinished = _progFile + "_finished";    
    ofstream ofs;
    
    // Create file if not yet there
    if (!boost::filesystem::exists(progFileFinished)) {
        ofs.open(progFileFinished.c_str(), ofstream::out);
        LOG(logERROR,*(thread->getLogger())) << "Created new file "
            << progFileFinished << flush;
    }
    else {
        ofs.open(progFileFinished.c_str(), ofstream::out | ofstream::app);
    }
    
    if (!ofs.is_open()) {
        LOG(logERROR,*(thread->getLogger())) << "Could not open file "
            << progFileFinished << flush;
        throw runtime_error("Bad file handle.");
    }
    
    JobItVec vjit;
    for (vjit = _jobsToSync.begin(); vjit != _jobsToSync.end(); ++vjit) {
        ofs << WriteProgLine(*vjit, "finished");
    }
    
    ofs.close();    
    _jobsToSync.clear();
    
    
    // READ PREVIOUSLY ASSIGNED JOBS INTO MAP
    string progFileAssigned = _progFile + "_assigned";
    
    // Create file if not yet there
    if (!boost::filesystem::exists(progFileAssigned)) {
        ofs.open(progFileAssigned.c_str(), ofstream::out);
        ofs.close();
        LOG(logERROR,*(thread->getLogger())) << "Created new file "
            << progFileAssigned << flush;
    }
    
    map<int, string> assignedId_line;

    ifstream ifs;
    ifs.open(progFileAssigned.c_str(), ifstream::in);

    if (!ifs.is_open()) {
        LOG(logERROR,*(thread->getLogger())) << "Could not open file "
            << progFileAssigned << flush;
        throw runtime_error("Bad file handle.");
    }    

    while ( ifs.good() ) {
        string line;
        std::getline(ifs, line);

        vector< string > split;
        Tokenizer toker(line, " \t");
        toker.ToVector(split);
        if ( !split.size()      ||
              split[0] == "#"   ||
              split[0].substr(0,1) == "#" ) { continue; }

        int id = boost::lexical_cast<int>(split[0]);
        assignedId_line[id] = line;
    }
    
    ifs.close();
    
    // GRAB NEW JOBS TO PROCESS
    _jobsToProc.clear();
    
    while(_jobsToProc.size() < _nThreads) {
        
        if (_metajit == _jobs->end()) break;
        
        
        if ( assignedId_line.count((*_metajit)->getId()) > 0 ) {
            ; // Already being processed somewhere
        }
        else {
            _jobsToProc.push_back(*_metajit);
        }
        
        ++_metajit;        
    }
    
    // APPEND NEWLY ASSIGNED JOBS TO LIST OF ASSIGNED JOBS
    ofs.open(progFileAssigned.c_str(), ofstream::out | ofstream::app);   
    
    if (!ofs.is_open()) {
        LOG(logERROR,*(thread->getLogger())) << "Could not open file "
            << progFileAssigned << flush;
        throw runtime_error("Bad file handle.");
    }
    
    for (vjit = _jobsToProc.begin(); vjit != _jobsToProc.end(); ++vjit) {
        ofs << WriteProgLine(*vjit, "assigned");
    }
    
    ofs.close();
    
    
    this->ReleaseProgFile(thread);
}


template<typename JobContainer, typename pJob>
void ProgObserver<JobContainer,pJob>::LockProgFile(QMThread *thread) {
    
    string sharedProgFolder = _progFile + "_locked";
    while (!boost::filesystem::create_directories(sharedProgFolder)) {
        LOG(logDEBUG,*(thread->getLogger())) 
            << "Waiting for lock on " << _progFile << flush;
        boost::this_thread::sleep(boost::posix_time::milliseconds(1000));
    }
    LOG(logDEBUG,*(thread->getLogger()))
        << "Imposed lock on " << _progFile << flush;
}


template<typename JobContainer, typename pJob>
void ProgObserver<JobContainer,pJob>::ReleaseProgFile(QMThread *thread) {
    
    string sharedProgFolder = _progFile + "_locked";
    LOG(logDEBUG,*(thread->getLogger()))
        << "Releasing " << _progFile << ". " << flush;
    boost::filesystem::remove(sharedProgFolder);
}


template<typename JobContainer, typename pJob>
string ProgObserver<JobContainer,pJob>::WriteProgLine(pJob job, string progStage) {
    
    string prog = "";
    if (progStage == "finished") prog = 
        (format("%1$5d %2$10s -FINISHED-") % job->getId() % job->getTag()).str();
    
    else if (progStage == "assigned") prog = 
        (format("%1$5d %2$10s -ASSIGNED-") % job->getId() % job->getTag()).str();
    
    else
        assert(false);
    
    return prog;    
}


template<typename JobContainer, typename pJob>
void ProgObserver<JobContainer,pJob>::ReportJobDone(pJob job, QMThread *thread) {
    
    _jobsToSync.push_back(job);
}


template<>
void ProgObserver< vector<XJob*>, XJob* >::ReportJobDone(XJob *job, QMThread *thread) {
    
    LOG(logDEBUG,*(thread->getLogger())) << "Reporting job done." << flush;
    _jobsToSync.push_back(job);
}


// REGISTER
template class ProgObserver< vector<XJob*>, XJob* >;
    
    
    
    
    
    
}}
