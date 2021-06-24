/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
/// For an earlier history see ctp repo commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

// Standard includes
#include <fstream>
#include <unistd.h>

// Third party includes
#include <boost/algorithm/string/replace.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <sys/types.h>

// Local VOTCA includes
#include "votca/xtp/job.h"
#include "votca/xtp/progressobserver.h"
#include "votca/xtp/qmthread.h"

using boost::format;

namespace votca {
namespace xtp {

template <typename JobContainer>
typename ProgObserver<JobContainer>::Job *
    ProgObserver<JobContainer>::RequestNextJob(QMThread &thread) {

  lockThread_.Lock();
  Job *jobToProc = nullptr;

  XTP_LOG(Log::error, thread.getLogger())
      << "Requesting next job" << std::flush;

  // NEED NEW CHUNK?
  if (nextjit_ == jobsToProc_.end() && moreJobsAvailable_) {
    SyncWithProgFile(thread);
    nextjit_ = jobsToProc_.begin();
    if (nextjit_ == jobsToProc_.end()) {
      moreJobsAvailable_ = false;
      XTP_LOG(Log::error, thread.getLogger())
          << "Sync did not yield any new jobs." << std::flush;
    }
  }

  // JOBS EATEN ALL UP?
  if (nextjit_ == jobsToProc_.end()) {
    if (maxJobs_ == startJobsCount_) {
      XTP_LOG(Log::error, thread.getLogger())
          << "Next job: ID = - (reached maximum for this process)"
          << std::flush;
    } else {
      XTP_LOG(Log::error, thread.getLogger())
          << "Next job: ID = - (none available)" << std::flush;
    }
  }
  // TAKE A BITE
  else {
    jobToProc = *nextjit_;
    ++nextjit_;
    XTP_LOG(Log::error, thread.getLogger())
        << "Next job: ID = " << jobToProc->getId() << std::flush;
  }

  if (!thread.isMaverick() && jobToProc != nullptr) {
    Index idx = jobToProc->getId();
    Index frac = (jobs_.size() >= 10) ? 10 : jobs_.size();
    Index rounded = Index(double(jobs_.size()) / double(frac)) * frac;
    Index tenth = rounded / frac;
    if (idx % tenth == 0) {
      double percent = double(idx) / double(rounded) * 100 + 0.5;
      std::cout << (format("=> [%1$2.0f%%] ") % percent).str() << std::flush;
    }
  }

  lockThread_.Unlock();
  return jobToProc;
}

template <typename JobContainer>
void ProgObserver<JobContainer>::ReportJobDone(Job &job, Result &res,
                                               QMThread &thread) {
  lockThread_.Lock();
  XTP_LOG(Log::error, thread.getLogger())
      << "Reporting job results" << std::flush;
  // RESULTS, TIME, HOST
  job.UpdateFromResult(res);
  job.setTime(GenerateTime());
  job.setHost(GenerateHost());
  // PRINT PROGRESS BAR
  jobsReported_ += 1;
  if (!thread.isMaverick()) {
    std::cout << std::endl << thread.getLogger() << std::flush;
  }
  lockThread_.Unlock();
  return;
}

template <typename JobContainer>
std::string ProgObserver<JobContainer>::GenerateHost() {
  char host[128];
  (void)gethostname(host, sizeof host);
  pid_t pid = getpid();
  return (format("%1$s:%2$d") % host % pid).str();
}

template <typename JobContainer>
std::string ProgObserver<JobContainer>::GenerateTime() {
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  return (format("%1$s") % now.time_of_day()).str();
}

template <typename JobContainer>
void ProgObserver<JobContainer>::SyncWithProgFile(QMThread &thread) {

  // INTERPROCESS FILE LOCKING (THREAD LOCK IN ::RequestNextJob)
  this->LockProgFile(thread);

  std::string progFile = progFile_;
  std::string progBackFile = progFile_ + "~";

  // LOAD EXTERNAL JOBS FROM SHARED XML & UPDATE INTERNAL JOBS
  XTP_LOG(Log::info, thread.getLogger())
      << "Update internal structures from job file" << std::flush;
  JobContainer jobs_ext = LOAD_JOBS(progFile);
  UPDATE_JOBS(jobs_ext, jobs_, GenerateHost());

  // GENERATE BACK-UP FOR SHARED XML
  XTP_LOG(Log::info, thread.getLogger())
      << "Create job-file back-up" << std::flush;
  WRITE_JOBS(jobs_, progBackFile);

  // ASSIGN NEW JOBS IF AVAILABLE
  XTP_LOG(Log::error, thread.getLogger())
      << "Assign jobs from stack" << std::flush;
  jobsToProc_.clear();

  Index cacheSize = cacheSize_;
  while (int(jobsToProc_.size()) < cacheSize) {
    if (metajit_ == jobs_.end() || startJobsCount_ == maxJobs_) {
      break;
    }

    bool startJob = false;

    // Start if job available or restart patterns matched
    if ((metajit_->isAvailable()) ||
        (restartMode_ && restart_stats_.count(metajit_->getStatusStr())) ||
        (restartMode_ && restart_hosts_.count(metajit_->getHost()))) {
      startJob = true;
    }

    if (startJob) {
      metajit_->Reset();
      metajit_->setStatus("ASSIGNED");
      metajit_->setHost(GenerateHost());
      metajit_->setTime(GenerateTime());
      jobsToProc_.push_back(&*metajit_);
      startJobsCount_ += 1;
    }

    ++metajit_;
  }

  // UPDATE PROGRESS STATUS FILE
  WRITE_JOBS(jobs_, progFile);

  // RELEASE PROGRESS STATUS FILE
  this->ReleaseProgFile(thread);
  return;
}

template <typename JobContainer>
void ProgObserver<JobContainer>::LockProgFile(QMThread &thread) {
  flock_ = std::unique_ptr<boost::interprocess::file_lock>(
      new boost::interprocess::file_lock(lockFile_.c_str()));
  flock_->lock();
  XTP_LOG(Log::warning, thread.getLogger())
      << "Imposed lock on " << lockFile_ << std::flush;
  XTP_LOG(Log::warning, thread.getLogger())
      << "Sleep ... " << lockFile_ << std::flush;
  XTP_LOG(Log::warning, thread.getLogger())
      << "Wake up ... " << lockFile_ << std::flush;
}

template <typename JobContainer>
void ProgObserver<JobContainer>::ReleaseProgFile(QMThread &thread) {
  flock_->unlock();
  XTP_LOG(Log::warning, thread.getLogger())
      << "Releasing " << lockFile_ << ". " << std::flush;
}

template <typename JobContainer>
void ProgObserver<JobContainer>::InitCmdLineOpts(
    const boost::program_options::variables_map &optsMap) {

  lockFile_ = optsMap["file"].as<std::string>();
  cacheSize_ = optsMap["cache"].as<Index>();
  maxJobs_ = optsMap["maxjobs"].as<Index>();
  std::string restartPattern = optsMap["restart"].as<std::string>();

  // restartPattern = e.g. host(pckr124:1234) stat(FAILED)
  boost::algorithm::replace_all(restartPattern, " ", "");
  if (restartPattern == "") {
    restartMode_ = false;
  } else {
    restartMode_ = true;
  }

  std::vector<std::string> patterns =
      tools::Tokenizer(restartPattern, "(,)").ToVector();

  std::string category = "";
  for (const std::string &pattern : patterns) {

    if (pattern == "host" || pattern == "stat") {
      category = pattern;

    } else if (category == "host") {
      restart_hosts_[pattern] = true;
    } else if (category == "stat") {
      if (pattern == "ASSIGNED" || pattern == "COMPLETE") {
        std::cout << "Restart if status == " << pattern
                  << "? Not necessarily a good idea." << std::endl;
      }
      restart_stats_[pattern] = true;
    }

    else {
      throw std::runtime_error(
          "Restart pattern ill-defined, format is"
          "[host([HOSTNAME:PID])] [stat([STATUS])]");
    }
  }
  return;
}

template <typename JobContainer>
void ProgObserver<JobContainer>::InitFromProgFile(std::string progFile,
                                                  QMThread &thread) {

  progFile_ = progFile;
  jobsReported_ = 0;

  XTP_LOG(Log::error, thread.getLogger())
      << "Job file = '" << progFile_ << "', ";
  XTP_LOG(Log::info, thread.getLogger())
      << "lock file = '" << lockFile_ << "', ";
  XTP_LOG(Log::error, thread.getLogger())
      << "cache size =  " << cacheSize_ << std::flush;

  XTP_LOG(Log::error, thread.getLogger())
      << "Initialize jobs from " << progFile << std::flush;
  XTP_LOG(Log::info, thread.getLogger()) << "Lock & load " << std::flush;

  // LOCK, READ INTO XML
  this->LockProgFile(thread);

  // ... Clear container
  jobs_.clear();

  // ... Load new, set availability bool
  jobs_ = LOAD_JOBS(progFile);
  metajit_ = jobs_.begin();
  WRITE_JOBS(jobs_, progFile + "~");
  XTP_LOG(Log::error, thread.getLogger())
      << "Registered " << jobs_.size() << " jobs." << std::flush;
  if (jobs_.size() > 0) {
    moreJobsAvailable_ = true;
  } else {
    moreJobsAvailable_ = false;
  }

  // SUMMARIZE OBSERVER VARIABLES: RESTART PATTERN, CACHE, LOCK FILE
  if (restartMode_ && restart_hosts_.size()) {
    std::string infostr = "Restart if host == ";
    for (const std::pair<const std::string, bool> &host : restart_hosts_) {
      infostr += host.first + " ";
    }
    XTP_LOG(Log::error, thread.getLogger()) << infostr << std::flush;
  }
  if (restartMode_ && restart_stats_.size()) {
    std::string infostr = "Restart if stat == ";
    for (const std::pair<const std::string, bool> &host : restart_hosts_) {
      infostr += host.first + " ";
    }
    XTP_LOG(Log::error, thread.getLogger()) << infostr << std::flush;
  }

  // RELEASE PROGRESS FILE
  this->ReleaseProgFile(thread);
  return;
}

// REGISTER
template class ProgObserver<std::vector<Job> >;

}  // namespace xtp
}  // namespace votca
