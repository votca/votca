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
#include <cstring>
#include <stdexcept>

// Local VOTCA includes
#include "votca/xtp/symmetric_tiles.h"

namespace votca {
namespace xtp {

double SymmetricTiles::Bytes(Index n, Index tile) {
  if (n <= 0) {
    return 0.0;
  }
  const Index nt = (n + tile - 1) / tile;
  const double last = double(n - (nt - 1) * tile);
  const double full = double(tile);
  // full tiles of the strict lower triangle of the first nt-1 tile rows,
  // the last tile row, and the diagonal tiles
  const double nfull_lower = double(nt - 1) * double(nt - 2) / 2.0;
  double elements = nfull_lower * full * full;
  elements += double(nt - 1) * last * full;                // last tile row
  elements += double(nt - 1) * full * full + last * last;  // diagonal tiles
  return 8.0 * elements;
}

void SymmetricTiles::Allocate(Index n, Index kind, Index tile) {
  Release();
  if (tile <= 0) {
    throw std::runtime_error("SymmetricTiles: tile size must be positive");
  }
  n_ = n;
  tile_ = tile;
  kind_ = kind;
  ntiles_ = (n + tile - 1) / tile;
  offset_.assign(std::size_t(ntiles_ * (ntiles_ + 1) / 2), 0);
  std::size_t total = 0;
  for (Index I = 0; I < ntiles_; ++I) {
    for (Index J = 0; J <= I; ++J) {
      offset_[std::size_t(I * (I + 1) / 2 + J)] = total;
      total += std::size_t(TileSize(I) * TileSize(J));
    }
  }
  // not value-initialised: the build writes every element (in parallel, so
  // the pages are spread over the threads that use them)
  data_.reset(new double[total]);
}

void SymmetricTiles::Release() {
  data_.reset();
  offset_.clear();
  n_ = 0;
  ntiles_ = 0;
  built_ = false;
}

void SymmetricTiles::SetRowLower(Index r, const double* values) {
  const Index I = r / tile_;
  const Index local = r - TileStart(I);
  for (Index J = 0; J <= I; ++J) {
    const Index j0 = TileStart(J);
    const Index len = (J < I) ? TileSize(J) : local + 1;
    double* dest = TileData(I, J) + local * TileSize(J);
    std::memcpy(dest, values + j0, std::size_t(len) * sizeof(double));
  }
}

void SymmetricTiles::FinishBuild() {
#pragma omp parallel for schedule(dynamic)
  for (Index I = 0; I < ntiles_; ++I) {
    RowMajorMap T(TileData(I, I), TileSize(I), TileSize(I));
    T.triangularView<Eigen::StrictlyUpper>() = T.transpose();
  }
  built_ = true;
}

Eigen::MatrixXd SymmetricTiles::Multiply(const Eigen::MatrixXd& X,
                                         double factor) const {
  if (!built_) {
    throw std::runtime_error("SymmetricTiles::Multiply: matrix not built");
  }
  if (X.rows() != n_) {
    throw std::runtime_error("SymmetricTiles::Multiply: size mismatch");
  }
  Eigen::MatrixXd Y(n_, X.cols());
  // each task owns one row block of Y: the tiles of tile row I and, for the
  // part above the diagonal, the transposes of the tiles of tile column I
#pragma omp parallel for schedule(dynamic)
  for (Index I = 0; I < ntiles_; ++I) {
    const Index i0 = TileStart(I);
    const Index bi = TileSize(I);
    auto Yi = Y.middleRows(i0, bi);
    Yi.setZero();
    for (Index J = 0; J <= I; ++J) {
      ConstRowMajorMap T(TileData(I, J), bi, TileSize(J));
      Yi.noalias() += T * X.middleRows(TileStart(J), TileSize(J));
    }
    for (Index J = I + 1; J < ntiles_; ++J) {
      ConstRowMajorMap T(TileData(J, I), TileSize(J), bi);
      Yi.noalias() += T.transpose() * X.middleRows(TileStart(J), TileSize(J));
    }
    if (factor != 1.0) {
      Yi *= factor;
    }
  }
  return Y;
}

double SymmetricTiles::operator()(Index i, Index j) const {
  if (i < j) {
    std::swap(i, j);
  }
  const Index I = i / tile_;
  const Index J = j / tile_;
  return TileData(I, J)[(i - TileStart(I)) * TileSize(J) + (j - TileStart(J))];
}

}  // namespace xtp
}  // namespace votca
