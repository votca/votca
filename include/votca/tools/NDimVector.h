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

#ifndef VOTCA_TOOLS_NDIMVECTOR_H
#define VOTCA_TOOLS_NDIMVECTOR_H

// Standard includes
#include "types.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
#include <vector>

namespace votca {
namespace tools {

/**
 * \brief N-Dim Vector
 *
 * This is basically a normal std::vector with an additional
 * operator(i,j,k,l...), which allows you to treat it as a higher dim object
 * It is column major (the first index is the fast one)
 *
 */

template <class T, int dim>
class NDimVector {

  static constexpr int dim_ = dim;

  template <typename... IndexTypes>
  Index linearIndex(Index firstDimension, IndexTypes... otherDimensions) const {
    static_assert((sizeof...(otherDimensions) + 1 == dim_),
                  "Number of dimensions given does not match rank");
    std::array<Index, dim_> indices = {firstDimension, otherDimensions...};
    assert(std::all_of(indices.begin(), indices.end(),
                       [](Index size) { return size >= 0; }) &&
           "All indeces must be non-negative");
    assert(std::equal(indices.begin(), indices.end(), dimensions_.begin(),
                      [](Index index, Index dimension) {
                        return index < dimension;
                      }) &&
           "All indeces must be smaller than dimension");
    return std::inner_product(indices.begin(), indices.end(), offsets_.begin(),
                              0);
  }

 public:
  constexpr Index rank() { return dimensions_.size(); }

  NDimVector() = default;

  template <typename... IndexTypes>
  NDimVector(Index firstDimension, IndexTypes... otherDimensions) {
    // The number of dimensions used to construct a tensor must be equal to the
    // rank of the tensor.
    static_assert((sizeof...(otherDimensions) + 1 == dim_),
                  "Number of dimensions given does not match rank");
    dimensions_ = {firstDimension, otherDimensions...};

    assert(std::all_of(dimensions_.begin(), dimensions_.end(),
                       [](Index size) { return size >= 0; }) &&
           "All dimensions must be larger 0");
    Index size = std::accumulate(dimensions_.begin(), dimensions_.end(), 1,
                                 std::multiplies<Index>());
    storage_.resize(size);
    offsets_[0] = 1;
    std::partial_sum(dimensions_.begin(), dimensions_.end() - 1,
                     offsets_.begin() + 1, std::multiplies<Index>());
  }

  template <typename... IndexTypes>
  T& operator()(Index firstDimension, IndexTypes... otherDimensions) {
    return storage_[linearIndex(firstDimension, otherDimensions...)];
  }

  template <typename... IndexTypes>
  const T& operator()(Index firstDimension,
                      IndexTypes... otherDimensions) const {
    return storage_[linearIndex(firstDimension, otherDimensions...)];
  }

  Index dimension(Index i) const { return dimensions_[i]; }
  Index size() const { return storage_.size(); }

  typename std::vector<T>::iterator begin() { return storage_.begin(); }
  typename std::vector<T>::iterator end() { return storage_.end(); }
  typename std::vector<T>::const_iterator begin() const {
    return storage_.begin();
  }
  typename std::vector<T>::const_iterator end() const { return storage_.end(); }

 private:
  std::vector<T> storage_;
  std::array<Index, dim> dimensions_;
  std::array<Index, dim> offsets_;
};
}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_NDIMVECTOR_H
