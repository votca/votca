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
#ifndef VOTCA_XTP_SCREENING_KERNELS_H
#define VOTCA_XTP_SCREENING_KERNELS_H

// Standard includes
#include <functional>
#include <optional>
#include <string>

// Local VOTCA includes
#include "dfttimings.h"
#include "eigen.h"

namespace votca {
namespace xtp {

/// Weighted Gram matrix S = sum_m R_m^T diag(w_m) R_m over a set of row
/// blocks R_m (n_m x ncols), as in the RPA response sum_{v,c} M_vc^T d_vc M_vc.
///
/// Only the lower triangle is computed (half the flops of the full product),
/// in tiles of S over row blocks stacked from several m, each tile one GEMM;
/// rows with negative weights are kept apart, so S = Xp^T Xp - Xn^T Xn with
/// X = sqrt|w| R. One output matrix, no per-thread copies of S.
class WeightedGram {
 public:
  /// Calls use(R_m); R_m may be a block of stored data or a temporary.
  using RowsVisitor =
      std::function<void(const Eigen::Ref<const Eigen::MatrixXd>&)>;
  using RowsFn = std::function<void(Index m, const RowsVisitor& use)>;
  /// Writes the weights of block m; their number is the number of rows.
  using WeightsFn = std::function<void(Index m, Eigen::VectorXd& w)>;

  /// Bytes for the stacked rows of one batch of blocks.
  static void setBatchMemory(double bytes) { batch_bytes_ = bytes; }

  static Eigen::MatrixXd Compute(Index nblocks, Index ncols,
                                 const WeightsFn& weights, const RowsFn& rows);

 private:
  static double batch_bytes_;
};

/// Eigenvalues (ascending) and eigenvectors of a symmetric matrix; uses the
/// divide-and-conquer LAPACK driver when built with MKL.
struct SymmetricEigenSystem {
  Eigen::VectorXd values;
  Eigen::MatrixXd vectors;
};
/// With vectors = false only the eigenvalues are computed.
SymmetricEigenSystem SymmetricEigen(const Eigen::MatrixXd& A,
                                    bool vectors = true);

/// Cholesky factor A = L L^T of a symmetric positive definite matrix, with
/// the same backends as SymmetricEigen (Accelerate LAPACK on macOS, LAPACKE
/// with MKL or Accelerate, otherwise Eigen). ok() is false if A is not
/// positive definite.
class CholeskyFactor {
 public:
  CholeskyFactor() = default;
  explicit CholeskyFactor(const Eigen::MatrixXd& A);
  void compute(const Eigen::MatrixXd& A);
  bool ok() const { return ok_; }
  const Eigen::MatrixXd& matrixL() const { return L_; }
  /// L^-1 B
  Eigen::MatrixXd SolveL(const Eigen::MatrixXd& B) const;
  /// A^-1 B
  Eigen::MatrixXd Solve(const Eigen::MatrixXd& B) const;

 private:
  Eigen::MatrixXd L_;
  bool ok_ = false;
};

/// Inverse of a symmetric positive definite matrix (Cholesky); falls back to
/// LU if the Cholesky factorisation fails.
Eigen::MatrixXd InverseSPD(const Eigen::MatrixXd& A);

/// A timing scope that is only opened when a DFTTimings object is given and
/// the caller is not inside an OpenMP parallel region (DFTTimings is not
/// thread safe; e.g. CDA builds epsilon inside the parallel QP solver).
class OptionalTiming {
 public:
  OptionalTiming(DFTTimings* timings, const std::string& name) {
    if (timings != nullptr && !OPENMP::InsideActiveParallelRegion()) {
      scope_.emplace(*timings, name);
    }
  }

 private:
  std::optional<DFTTimings::Scope> scope_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SCREENING_KERNELS_H
