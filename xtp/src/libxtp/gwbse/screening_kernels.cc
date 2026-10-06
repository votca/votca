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
#include <cmath>
#include <utility>
#include <vector>

// Local VOTCA includes
#include "votca/xtp/screening_kernels.h"

namespace votca {
namespace xtp {

double WeightedGram::batch_bytes_ = 5e8;

namespace {

// Lower-triangle tiles (I >= J) of an n x n matrix in blocks of size bs.
std::vector<std::pair<Index, Index>> LowerTiles(Index n, Index bs) {
  const Index nt = (n + bs - 1) / bs;
  std::vector<std::pair<Index, Index>> tiles;
  tiles.reserve(std::size_t(nt * (nt + 1) / 2));
  for (Index J = 0; J < nt; ++J) {
    for (Index I = J; I < nt; ++I) {
      tiles.emplace_back(I, J);
    }
  }
  return tiles;
}

// S_lower += Xp^T Xp - Xn^T Xn on the given tiles, one tile per task.
void AccumulateTiles(Eigen::MatrixXd& S, const Eigen::MatrixXd& Xp,
                     const Eigen::MatrixXd& Xn, Index bs,
                     const std::vector<std::pair<Index, Index>>& tiles) {
  const Index n = S.cols();
#pragma omp parallel for schedule(dynamic)
  for (std::size_t t = 0; t < tiles.size(); ++t) {
    const Index i0 = tiles[t].first * bs;
    const Index j0 = tiles[t].second * bs;
    const Index ni = std::min(bs, n - i0);
    const Index nj = std::min(bs, n - j0);
    auto block = S.block(i0, j0, ni, nj);
    if (Xp.rows() > 0) {
      block.noalias() +=
          Xp.middleCols(i0, ni).transpose() * Xp.middleCols(j0, nj);
    }
    if (Xn.rows() > 0) {
      block.noalias() -=
          Xn.middleCols(i0, ni).transpose() * Xn.middleCols(j0, nj);
    }
  }
}

}  // namespace

Eigen::MatrixXd WeightedGram::Compute(Index nblocks, Index ncols,
                                      const WeightsFn& weights,
                                      const RowsFn& rows) {
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(ncols, ncols);
  if (nblocks == 0 || ncols == 0) {
    return S;
  }
  // Tiles small enough to give every thread several tasks, large enough
  // for efficient GEMMs.
  const Index nthreads = OPENMP::getMaxThreads();
  Index bs = 256;
  while (bs > 64) {
    const Index nt = (ncols + bs - 1) / bs;
    if (nt * (nt + 1) / 2 >= 4 * nthreads) {
      break;
    }
    bs /= 2;
  }
  const auto tiles = LowerTiles(ncols, bs);

  // Called from inside a parallel region (e.g. per CDA residue), every
  // thread builds its own matrix: keep the stacked rows small then.
  const double budget = OPENMP::InsideActiveParallelRegion()
                            ? std::min(batch_bytes_, 6.4e7)
                            : batch_bytes_;
  const Index max_rows =
      std::max<Index>(1, static_cast<Index>(budget / (8.0 * double(ncols))));

  std::vector<Eigen::VectorXd> w(static_cast<std::size_t>(nblocks));
  Index first = 0;
  while (first < nblocks) {
    // collect blocks until the stacked rows reach the budget
    Index last = first;
    Index nrows = 0;
    do {
      weights(last, w[std::size_t(last)]);
      nrows += w[std::size_t(last)].size();
      ++last;
    } while (last < nblocks && nrows < max_rows);

    // rows with positive and with negative weight, per block
    std::vector<Index> pos_offset(std::size_t(last - first) + 1, 0);
    std::vector<Index> neg_offset(std::size_t(last - first) + 1, 0);
    for (Index m = first; m < last; ++m) {
      const Eigen::VectorXd& wm = w[std::size_t(m)];
      const Index npos = (wm.array() > 0.0).count();
      const Index nneg = (wm.array() < 0.0).count();
      const std::size_t k = std::size_t(m - first);
      pos_offset[k + 1] = pos_offset[k] + npos;
      neg_offset[k + 1] = neg_offset[k] + nneg;
    }
    Eigen::MatrixXd Xp(pos_offset.back(), ncols);
    Eigen::MatrixXd Xn(neg_offset.back(), ncols);

#pragma omp parallel for schedule(dynamic)
    for (Index m = first; m < last; ++m) {
      const std::size_t k = std::size_t(m - first);
      const Eigen::VectorXd& wm = w[std::size_t(m)];
      const Index npos = pos_offset[k + 1] - pos_offset[k];
      const Index nneg = neg_offset[k + 1] - neg_offset[k];
      auto Xp_m = Xp.middleRows(pos_offset[k], npos);
      auto Xn_m = Xn.middleRows(neg_offset[k], nneg);
      if (npos == wm.size()) {
        // all weights positive (epsilon at zero or imaginary frequency)
        const Eigen::VectorXd scale = wm.cwiseSqrt();
        rows(m, [&](const Eigen::Ref<const Eigen::MatrixXd>& R) {
          Xp_m.noalias() = scale.asDiagonal() * R;
        });
        continue;
      }
      std::vector<Index> pos_rows;
      std::vector<Index> neg_rows;
      pos_rows.reserve(std::size_t(npos));
      neg_rows.reserve(std::size_t(nneg));
      for (Index r = 0; r < wm.size(); ++r) {
        if (wm(r) > 0.0) {
          pos_rows.push_back(r);
        } else if (wm(r) < 0.0) {
          neg_rows.push_back(r);
        }
      }
      Eigen::VectorXd pos_scale(npos);
      Eigen::VectorXd neg_scale(nneg);
      for (Index i = 0; i < npos; ++i) {
        pos_scale(i) = std::sqrt(wm(pos_rows[std::size_t(i)]));
      }
      for (Index i = 0; i < nneg; ++i) {
        neg_scale(i) = std::sqrt(-wm(neg_rows[std::size_t(i)]));
      }
      rows(m, [&](const Eigen::Ref<const Eigen::MatrixXd>& R) {
        for (Index c = 0; c < ncols; ++c) {
          for (Index i = 0; i < npos; ++i) {
            Xp_m(i, c) = pos_scale(i) * R(pos_rows[std::size_t(i)], c);
          }
          for (Index i = 0; i < nneg; ++i) {
            Xn_m(i, c) = neg_scale(i) * R(neg_rows[std::size_t(i)], c);
          }
        }
      });
    }
#if defined(EIGEN_USE_MKL_ALL)
    if (!OPENMP::InsideActiveParallelRegion()) {
      // threaded MKL dsyrk on the lower triangle
      if (Xp.rows() > 0) {
        S.selfadjointView<Eigen::Lower>().rankUpdate(Xp.transpose(), 1.0);
      }
      if (Xn.rows() > 0) {
        S.selfadjointView<Eigen::Lower>().rankUpdate(Xn.transpose(), -1.0);
      }
      first = last;
      continue;
    }
#endif
    AccumulateTiles(S, Xp, Xn, bs, tiles);
    first = last;
  }
  S.triangularView<Eigen::StrictlyUpper>() = S.transpose();
  return S;
}

SymmetricEigenSystem SymmetricEigen(const Eigen::MatrixXd& A) {
  SymmetricEigenSystem result;
#if defined(EIGEN_USE_MKL_ALL) || defined(EIGEN_USE_LAPACKE)
  const lapack_int n = static_cast<lapack_int>(A.rows());
  result.vectors = A;
  result.values.resize(A.rows());
  if (n > 0) {
    const lapack_int info =
        LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'L', n, result.vectors.data(),
                       static_cast<lapack_int>(result.vectors.outerStride()),
                       result.values.data());
    if (info == 0) {
      return result;
    }
  }
#endif
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  result.values = es.eigenvalues();
  result.vectors = es.eigenvectors();
  return result;
}

Eigen::MatrixXd InverseSPD(const Eigen::MatrixXd& A) {
#if defined(EIGEN_USE_MKL_ALL) || defined(EIGEN_USE_LAPACKE)
  const lapack_int n = static_cast<lapack_int>(A.rows());
  Eigen::MatrixXd inv = A;
  if (n == 0) {
    return inv;
  }
  const lapack_int lda = static_cast<lapack_int>(inv.outerStride());
  if (LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', n, inv.data(), lda) == 0 &&
      LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', n, inv.data(), lda) == 0) {
    inv.triangularView<Eigen::StrictlyUpper>() = inv.transpose();
    return inv;
  }
#else
  Eigen::LLT<Eigen::MatrixXd> llt(A);
  if (llt.info() == Eigen::Success) {
    return llt.solve(Eigen::MatrixXd::Identity(A.rows(), A.cols()));
  }
#endif
  return A.partialPivLu().inverse();
}

}  // namespace xtp
}  // namespace votca
