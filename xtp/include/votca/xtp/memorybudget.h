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

#pragma once
#ifndef VOTCA_XTP_MEMORYBUDGET_H
#define VOTCA_XTP_MEMORYBUDGET_H

// Standard includes
#include <string>

// VOTCA includes
#include <votca/tools/types.h>

namespace votca {
namespace xtp {

/**
 * \brief Memory the user grants the process (--memory, GB), for choices that
 * trade memory for speed (caches, batch sizes).
 *
 * The budget is the total for the process; with several jobs running
 * concurrently in one process (xtp_parallel -t) each gets an equal share.
 * It is advisory, as %mem in Gaussian or %maxcore in ORCA: it steers optional
 * allocations, it does not limit the memory a calculation needs anyway.
 *
 * Without a budget a fixed allowance of DefaultAllowanceBytes() is available
 * for optional buffers, regardless of what the process already holds.
 */
class MemoryBudget {
 public:
  /// total_gb <= 0: no budget given
  static void Set(double total_gb, Index concurrent_jobs = 1);
  static void SetConcurrentJobs(Index concurrent_jobs);
  static bool Given() { return total_bytes_ > 0.0; }
  static double DefaultAllowanceBytes() { return 4e9; }

  /// The share of the budget of one job (0 without a budget)
  static double JobBudgetBytes();

  /// Bytes one job may still allocate for optional buffers: its share of the
  /// budget minus its share of the resident memory of the process, minus a
  /// margin (5% of the share, at least 1 GB); without a budget the default
  /// allowance.
  static double AvailableBytes();

  /// Resident memory of the process (-1 if unknown on this platform)
  static double ResidentBytes();

  /// e.g. "budget 400 GB, 46.1 GB in use" or "no --memory given"
  static std::string Describe();

 private:
  static double total_bytes_;
  static Index concurrent_jobs_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MEMORYBUDGET_H
