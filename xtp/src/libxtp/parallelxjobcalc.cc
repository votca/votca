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

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "xtp_libint2.h"

// Local VOTCA includes
#include "votca/xtp/parallelxjobcalc.h"

using boost::format;

namespace votca {
namespace xtp {

template <typename JobContainer>
bool ParallelXJobCalc<JobContainer>::Evaluate(const Topology &top) {
  libint2::initialize();
  // INITIALIZE PROGRESS OBSERVER
  std::string progFile = jobfile_;
  std::unique_ptr<JobOperator> master = std::unique_ptr<JobOperator>(
      new JobOperator(-1, top, *this, openmp_threads_));
  master->getLogger().setReportLevel(Log::current_level);
  master->getLogger().setMultithreading(true);
  master->getLogger().setPreface(Log::info, "\nMST INF");
  master->getLogger().setPreface(Log::error, "\nMST ERR");
  master->getLogger().setPreface(Log::warning, "\nMST WAR");
  master->getLogger().setPreface(Log::debug, "\nMST DBG");
  progObs_->InitFromProgFile(progFile, *(master.get()));

  // CREATE + EXECUTE THREADS (XJOB HANDLERS)
  std::vector<std::unique_ptr<JobOperator>> jobOps;

  for (Index id = 0; id < nThreads_; id++) {
    jobOps.push_back(std::unique_ptr<JobOperator>(
        new JobOperator(id, top, *this, openmp_threads_)));
  }

  for (Index id = 0; id < nThreads_; ++id) {
    CustomizeLogger(*jobOps[id]);
  }

  if (!maverick_) {
    std::cout << std::endl;  // REQUIRED FOR PROGRESS BAR IN OBSERVER
  }

  for (Index id = 0; id < nThreads_; id++) {
    jobOps[id]->Start();
  }

  for (Index id = 0; id < nThreads_; id++) {
    jobOps[id]->WaitDone();
  }

  if (!maverick_) {
    for (Index id = 0; id < nThreads_; id++) {
      std::cout << std::endl << (jobOps[id]->getLogger()) << std::flush;
    }
  }

  jobOps.clear();

  // SYNC REMAINING COMPLETE JOBS
  progObs_->SyncWithProgFile(*(master.get()));
  libint2::finalize();
  return true;
}

template <typename JobContainer>
void ParallelXJobCalc<JobContainer>::JobOperator::Run() {
  OPENMP::setMaxThreads(openmp_threads_);
  while (true) {
    Job *job = master_.progObs_->RequestNextJob(*this);

    if (job == nullptr) {
      break;
    } else {
      Result res = this->master_.EvalJob(top_, *job, *this);
      this->master_.progObs_->ReportJobDone(*job, res, *this);
    }
  }
}

template <typename JobContainer>
void ParallelXJobCalc<JobContainer>::ParseCommonOptions(
    const tools::Property &options) {
  std::cout << "\n... ... Initialized with " << nThreads_ << " threads.\n";

  maverick_ = (nThreads_ == 1) ? true : false;

  std::cout << "\n... ... Using " << openmp_threads_ << " openmp threads for "
            << nThreads_ << "x" << openmp_threads_ << "="
            << nThreads_ * openmp_threads_ << " total threads." << std::flush;
  jobfile_ = options.get(".job_file").as<std::string>();
  mapfile_ = options.get(".map_file").as<std::string>();
}

template <typename JobContainer>
void ParallelXJobCalc<JobContainer>::CustomizeLogger(QMThread &thread) {

  // CONFIGURE LOGGER
  Logger &log = thread.getLogger();
  log.setReportLevel(Log::current_level);
  log.setMultithreading(maverick_);

  log.setPreface(Log::info,
                 (format("\nT%1$02d INF ...") % thread.getId()).str());
  log.setPreface(Log::error,
                 (format("\nT%1$02d ERR ...") % thread.getId()).str());
  log.setPreface(Log::warning,
                 (format("\nT%1$02d WAR ...") % thread.getId()).str());
  log.setPreface(Log::debug,
                 (format("\nT%1$02d DBG ...") % thread.getId()).str());
}

// REGISTER PARALLEL CALCULATORS
template class ParallelXJobCalc<std::vector<Job>>;

}  // namespace xtp
}  // namespace votca
