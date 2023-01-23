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

#pragma once
#ifndef VOTCA_XTP_PROGRESSOBSERVER_H
#define VOTCA_XTP_PROGRESSOBSERVER_H

// Standard includes
#include <vector>

// Third party includes
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/program_options.hpp>

// VOTCA includes
#include <votca/tools/mutex.h>
#include <votca/tools/property.h>

namespace votca {
namespace xtp {

class QMThread;

template <typename JobContainer>
class ProgObserver {

  using Job = typename JobContainer::value_type;
  using Result = typename Job::JobResult;

 public:
  void InitCmdLineOpts(const boost::program_options::variables_map &optsMap);
  void InitFromProgFile(std::string progFile, QMThread &thread);
  ProgObserver::Job *RequestNextJob(QMThread &thread);
  void ReportJobDone(Job &job, Result &res, QMThread &thread);

  void SyncWithProgFile(QMThread &thread);
  void LockProgFile(QMThread &thread);
  void ReleaseProgFile(QMThread &thread);

  std::string GenerateHost();
  std::string GenerateTime();

 private:
  std::string lockFile_ = "";
  std::string progFile_ = "";
  Index cacheSize_ = -1;
  JobContainer jobs_;

  std::vector<Job *> jobsToProc_;
  std::vector<Job *> jobsToSync_;

  using iterator = typename JobContainer::iterator;
  iterator metajit_;
  using iterator_vec = typename std::vector<Job *>::iterator;
  iterator_vec nextjit_;
  tools::Mutex lockThread_;
  std::unique_ptr<boost::interprocess::file_lock> flock_;

  std::map<std::string, bool> restart_hosts_;
  std::map<std::string, bool> restart_stats_;
  bool restartMode_ = false;
  Index jobsReported_ = 0;

  bool moreJobsAvailable_ = false;
  Index startJobsCount_ = 0;
  Index maxJobs_ = 0;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_PROGRESSOBSERVER_H
