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
#ifndef VOTCA_XTP_GRIDBOXACCUMULATOR_H
#define VOTCA_XTP_GRIDBOXACCUMULATOR_H

// Standard includes
#include <cstddef>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// Local VOTCA includes
#include "votca/xtp/eigen.h"
#include "votca/xtp/gridbox.h"

namespace votca {
namespace xtp {

/**
 * \brief Adds per-grid-box blocks into a set of AO matrices from inside an
 * OpenMP loop over grid boxes.
 *
 * A single critical section around all matrices serializes the loop when
 * there are many of them (one per response density or per atom). Two
 * strategies instead:
 *  - thread-local copies of all target matrices, summed once in Finish(), when
 *    they fit into the memory budget (no synchronization at all);
 *  - otherwise one lock per matrix, taken in a thread-dependent order so that
 *    threads work on different matrices.
 * Results are identical up to the order of floating-point additions.
 *
 * Usage: construct before the parallel loop, call Add() inside it (any
 * thread), call Finish() after it.
 */
class GridBoxAccumulator {
 public:
  /// targets: zero-initialized matrices of equal size; budget in bytes for
  /// the thread-local copies (all threads together)
  explicit GridBoxAccumulator(std::vector<Eigen::MatrixXd*> targets,
                              double budget_bytes = 2.0e9)
      : targets_(std::move(targets)) {
    const Index n = Index(targets_.size());
    const Index nthreads = OPENMP::getMaxThreads();
    double bytes = 0.0;
    for (const Eigen::MatrixXd* m : targets_) {
      bytes += double(m->size()) * double(sizeof(double));
    }
    thread_local_ =
        nthreads > 1 && double(nthreads - 1) * bytes <= budget_bytes;
    if (thread_local_) {
      // thread 0 writes into the targets themselves; the others get copies on
      // first use
      copies_.resize(std::size_t(nthreads));
    }
#ifdef _OPENMP
    if (!thread_local_ && nthreads > 1) {
      locks_.resize(std::size_t(n));
      for (omp_lock_t& l : locks_) {
        omp_init_lock(&l);
      }
    }
#endif
    (void)n;
  }

  ~GridBoxAccumulator() {
#ifdef _OPENMP
    for (omp_lock_t& l : locks_) {
      omp_destroy_lock(&l);
    }
#endif
  }
  GridBoxAccumulator(const GridBoxAccumulator&) = delete;
  GridBoxAccumulator& operator=(const GridBoxAccumulator&) = delete;

  bool UsesThreadLocalCopies() const { return thread_local_; }

  /// adds the box blocks small[k] into target k, for all k (same order as the
  /// targets); entries with size 0 are skipped
  void Add(const GridBox& box,
           const std::vector<const Eigen::MatrixXd*>& small) {
    const Index n = Index(targets_.size());
    const Index tid = OPENMP::getThreadId();
    if (thread_local_) {
      if (tid == 0) {
        for (Index k = 0; k < n; ++k) {
          if (small[std::size_t(k)]->size() > 0) {
            box.AddtoBigMatrix(*targets_[std::size_t(k)],
                               *small[std::size_t(k)]);
          }
        }
        return;
      }
      if (tid >= Index(copies_.size())) {  // more threads than at construction
#pragma omp critical(votca_xtp_gridboxaccumulator)
        {
          for (Index k = 0; k < n; ++k) {
            if (small[std::size_t(k)]->size() > 0) {
              box.AddtoBigMatrix(*targets_[std::size_t(k)],
                                 *small[std::size_t(k)]);
            }
          }
        }
        return;
      }
      std::vector<Eigen::MatrixXd>& mine = copies_[std::size_t(tid)];
      if (mine.empty()) {
        mine.reserve(std::size_t(n));
        for (const Eigen::MatrixXd* m : targets_) {
          mine.push_back(Eigen::MatrixXd::Zero(m->rows(), m->cols()));
        }
      }
      for (Index k = 0; k < n; ++k) {
        if (small[std::size_t(k)]->size() > 0) {
          box.AddtoBigMatrix(mine[std::size_t(k)], *small[std::size_t(k)]);
        }
      }
      return;
    }
#ifdef _OPENMP
    if (!locks_.empty()) {
      for (Index j = 0; j < n; ++j) {
        const Index k = (j + tid) % n;
        if (small[std::size_t(k)]->size() == 0) {
          continue;
        }
        omp_set_lock(&locks_[std::size_t(k)]);
        box.AddtoBigMatrix(*targets_[std::size_t(k)], *small[std::size_t(k)]);
        omp_unset_lock(&locks_[std::size_t(k)]);
      }
      return;
    }
#endif
    for (Index k = 0; k < n; ++k) {
      if (small[std::size_t(k)]->size() > 0) {
        box.AddtoBigMatrix(*targets_[std::size_t(k)], *small[std::size_t(k)]);
      }
    }
  }

  /// sums the thread-local copies into the targets (call outside the loop)
  void Finish() {
    if (!thread_local_) {
      return;
    }
    const Index n = Index(targets_.size());
#pragma omp parallel for schedule(static)
    for (Index k = 0; k < n; ++k) {
      for (const std::vector<Eigen::MatrixXd>& c : copies_) {
        if (!c.empty()) {
          *targets_[std::size_t(k)] += c[std::size_t(k)];
        }
      }
    }
    copies_.clear();
  }

 private:
  std::vector<Eigen::MatrixXd*> targets_;
  bool thread_local_ = false;
  std::vector<std::vector<Eigen::MatrixXd>> copies_;
#ifdef _OPENMP
  std::vector<omp_lock_t> locks_;
#endif
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GRIDBOXACCUMULATOR_H
