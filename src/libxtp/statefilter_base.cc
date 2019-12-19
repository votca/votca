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
#include <votca/xtp/statefilter_base.h>
namespace votca {
namespace xtp {

template <bool larger>
std::vector<Index> StateFilter_base::ReduceAndSortIndeces(
    const Eigen::VectorXd& overlap, Index offset, double threshold) const {
  std::vector<Index> index = std::vector<Index>(overlap.size());
  std::iota(index.begin(), index.end(), offset);
  std::stable_sort(index.begin(), index.end(), [&overlap](Index i1, Index i2) {
    if (larger) {
      return overlap[i1] > overlap[i2];
    } else {
      return overlap[i1] < overlap[i2];
    }
  });
  Index validelements;
  if (larger) {
    validelements = (overlap.array() > threshold).sum();
  } else {
    validelements = (overlap.array() < threshold).sum();
  }

  std::vector<Index> indexes(validelements);
  std::copy_n(index.begin(), validelements, indexes.begin());
  return indexes;
}

std::vector<Index> StateFilter_base::ReduceAndSortIndecesUp(
    const Eigen::VectorXd& overlap, Index offset, double threshold) const {
  return ReduceAndSortIndeces<true>(overlap, offset, threshold);
}
std::vector<Index> StateFilter_base::ReduceAndSortIndecesDown(
    const Eigen::VectorXd& overlap, Index offset, double threshold) const {
  return ReduceAndSortIndeces<false>(overlap, offset, threshold);
}

}  // namespace xtp
}  // namespace votca