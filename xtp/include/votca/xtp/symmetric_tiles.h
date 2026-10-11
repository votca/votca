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
#ifndef VOTCA_XTP_SYMMETRIC_TILES_H
#define VOTCA_XTP_SYMMETRIC_TILES_H

// Standard includes
#include <algorithm>
#include <memory>
#include <vector>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief A real symmetric n x n matrix stored as the tiles of its lower
 * triangle (about half the memory of the dense matrix).
 *
 * The tiles (I, J), I >= J, are tile x tile blocks (smaller at the end),
 * row-major, so that a row segment within a tile is contiguous. The diagonal
 * tiles are stored in full. Rows are written with SetRowLower, which takes
 * the columns 0..r of row r; FinishBuild() then fills the upper halves of the
 * diagonal tiles.
 *
 * Multiply(X) gives S X with the tiles of the lower triangle and their
 * transposes, one row block of the result per task.
 */
class SymmetricTiles {
 public:
  static constexpr Index DefaultTile = 1024;

  SymmetricTiles() = default;

  /// Bytes the tiles of an n x n matrix need
  static double Bytes(Index n, Index tile = DefaultTile);

  /// Allocates (uninitialised) storage for an n x n matrix; kind is a tag of
  /// what the matrix holds, checked by users that share it.
  void Allocate(Index n, Index kind, Index tile = DefaultTile);
  void Release();

  bool allocated() const { return data_ != nullptr; }
  bool built() const { return built_; }
  Index size() const { return n_; }
  Index kind() const { return kind_; }

  /// Row r, columns 0..r (values[0..r]); rows may be set in parallel.
  void SetRowLower(Index r, const double* values);
  /// Fills the upper halves of the diagonal tiles; marks the matrix built.
  void FinishBuild();

  /// factor * S * X
  Eigen::MatrixXd Multiply(const Eigen::MatrixXd& X, double factor = 1.0) const;

  /// Element (i, j), for tests
  double operator()(Index i, Index j) const;

 private:
  using RowMajorMap = Eigen::Map<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;
  using ConstRowMajorMap =
      Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>>;

  Index TileStart(Index I) const { return I * tile_; }
  Index TileSize(Index I) const { return std::min(tile_, n_ - TileStart(I)); }
  double* TileData(Index I, Index J) const {
    return data_.get() + offset_[std::size_t(I * (I + 1) / 2 + J)];
  }

  Index n_ = 0;
  Index tile_ = DefaultTile;
  Index ntiles_ = 0;
  Index kind_ = 0;
  bool built_ = false;
  std::vector<std::size_t> offset_;
  std::unique_ptr<double[]> data_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SYMMETRIC_TILES_H
