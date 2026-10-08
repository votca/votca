/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

// Standard includes
#include <algorithm>
#include <fstream>
#include <sstream>

#if defined(__APPLE__)
#include <mach/mach.h>
#endif

// Local VOTCA includes
#include "votca/xtp/memorybudget.h"

namespace votca {
namespace xtp {

double MemoryBudget::total_bytes_ = 0.0;
Index MemoryBudget::concurrent_jobs_ = 1;

void MemoryBudget::Set(double total_gb, Index concurrent_jobs) {
  total_bytes_ = (total_gb > 0.0) ? total_gb * 1e9 : 0.0;
  SetConcurrentJobs(concurrent_jobs);
}

void MemoryBudget::SetConcurrentJobs(Index concurrent_jobs) {
  concurrent_jobs_ = std::max<Index>(1, concurrent_jobs);
}

double MemoryBudget::JobBudgetBytes() {
  return total_bytes_ / double(concurrent_jobs_);
}

double MemoryBudget::ResidentBytes() {
#if defined(__APPLE__)
  mach_task_basic_info_data_t info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS) {
    return double(info.resident_size);
  }
  return -1.0;
#else
  std::ifstream status("/proc/self/status");
  std::string word;
  while (status >> word) {
    if (word == "VmRSS:") {
      double kb = 0.0;
      status >> kb;
      return kb * 1024.0;
    }
  }
  return -1.0;
#endif
}

double MemoryBudget::AvailableBytes() {
  if (!Given()) {
    return DefaultAllowanceBytes();
  }
  const double share = JobBudgetBytes();
  const double resident = std::max(0.0, ResidentBytes());
  const double margin = std::max(1e9, 0.05 * share);
  return std::max(0.0, share - resident / double(concurrent_jobs_) - margin);
}

std::string MemoryBudget::Describe() {
  std::ostringstream out;
  out.precision(3);
  if (!Given()) {
    out << "no --memory given, " << DefaultAllowanceBytes() * 1e-9
        << " GB for optional buffers";
    return out.str();
  }
  out << "budget " << JobBudgetBytes() * 1e-9 << " GB";
  if (concurrent_jobs_ > 1) {
    out << " per job (" << total_bytes_ * 1e-9 << " GB for " << concurrent_jobs_
        << " jobs)";
  }
  const double resident = ResidentBytes();
  if (resident >= 0.0) {
    out << ", " << resident * 1e-9 << " GB in use";
  }
  return out.str();
}

}  // namespace xtp
}  // namespace votca
