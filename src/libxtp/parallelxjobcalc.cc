/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <votca/xtp/parallelxjobcalc.h>

using boost::format;

namespace votca {
namespace xtp {

template <typename JobContainer>
bool ParallelXJobCalc<JobContainer>::EvaluateFrame(const Topology &top) {

  // INITIALIZE PROGRESS OBSERVER
  std::string progFile = _jobfile;
  std::unique_ptr<JobOperator> master = std::unique_ptr<JobOperator>(
      new JobOperator(-1, top, *this, _openmp_threads));
  master->getLogger().setReportLevel(logDEBUG);
  master->getLogger().setMultithreading(true);
  master->getLogger().setPreface(logINFO, "\nMST INF");
  master->getLogger().setPreface(logERROR, "\nMST ERR");
  master->getLogger().setPreface(logWARNING, "\nMST WAR");
  master->getLogger().setPreface(logDEBUG, "\nMST DBG");
  _progObs->InitFromProgFile(progFile, *(master.get()));

  // CREATE + EXECUTE THREADS (XJOB HANDLERS)
  std::vector<std::unique_ptr<JobOperator>> jobOps;

  for (unsigned int id = 0; id < _nThreads; id++) {
    jobOps.push_back(std::unique_ptr<JobOperator>(
        new JobOperator(id, top, *this, _openmp_threads)));
  }

  for (unsigned int id = 0; id < _nThreads; ++id) {
    CustomizeLogger(*jobOps[id]);
  }

  if (!_maverick) {
    std::cout << std::endl;  // REQUIRED FOR PROGRESS BAR IN OBSERVER
  }

  for (unsigned int id = 0; id < _nThreads; id++) {
    jobOps[id]->Start();
  }

  for (unsigned int id = 0; id < _nThreads; id++) {
    jobOps[id]->WaitDone();
  }

  if (!_maverick) {
    for (unsigned int id = 0; id < _nThreads; id++) {
      std::cout << std::endl << (jobOps[id]->getLogger()) << std::flush;
    }
  }

  jobOps.clear();

  // SYNC REMAINING COMPLETE JOBS
  _progObs->SyncWithProgFile(*(master.get()));

  return true;
}

template <typename JobContainer>
void ParallelXJobCalc<JobContainer>::JobOperator::Run() {
  OPENMP::setMaxThreads(_openmp_threads);
  while (true) {
    Job *job = _master._progObs->RequestNextJob(*this);

    if (job == nullptr) {
      break;
    } else {
      Result res = this->_master.EvalJob(_top, *job, *this);
      this->_master._progObs->ReportJobDone(*job, res, *this);
    }
  }
}

template <typename JobContainer>
void ParallelXJobCalc<JobContainer>::ParseCommonOptions(
    const tools::Property &options) {
  std::cout << std::endl
            << "... ... Initialized with " << _nThreads << " threads. "
            << std::flush;

  _maverick = (_nThreads == 1) ? true : false;

  std::string key = "options." + Identify();
  _openmp_threads = options.ifExistsReturnElseReturnDefault<int>(
      key + ".openmp_threads", _openmp_threads);
  std::cout << std::endl
            << "... ... Using " << _openmp_threads << " openmp threads for "
            << _nThreads << "x" << _openmp_threads << "="
            << _nThreads * _openmp_threads << " total threads." << std::flush;
  OPENMP::setMaxThreads(_openmp_threads);
  _jobfile = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".job_file");
  _mapfile = options.ifExistsReturnElseThrowRuntimeError<std::string>(
      key + ".map_file");
}

template <typename JobContainer>
void ParallelXJobCalc<JobContainer>::CustomizeLogger(QMThread &thread) {

  // CONFIGURE LOGGER
  Logger &log = thread.getLogger();
  log.setReportLevel(logDEBUG);
  log.setMultithreading(_maverick);

  log.setPreface(logINFO, (format("\nT%1$02d INF ...") % thread.getId()).str());
  log.setPreface(logERROR,
                 (format("\nT%1$02d ERR ...") % thread.getId()).str());
  log.setPreface(logWARNING,
                 (format("\nT%1$02d WAR ...") % thread.getId()).str());
  log.setPreface(logDEBUG,
                 (format("\nT%1$02d DBG ...") % thread.getId()).str());
}

// REGISTER PARALLEL CALCULATORS
template class ParallelXJobCalc<std::vector<Job>>;

}  // namespace xtp
}  // namespace votca
