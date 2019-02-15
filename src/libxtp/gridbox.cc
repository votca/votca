/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <votca/xtp/gridbox.h>

namespace votca {
namespace xtp {

void GridBox::AddtoBigMatrix(Eigen::MatrixXd& bigmatrix,
                             const Eigen::MatrixXd& smallmatrix) const {
  for (unsigned i = 0; i < ranges.size(); i++) {
    for (unsigned j = 0; j < ranges.size(); j++) {
      bigmatrix.block(ranges[i].start, ranges[j].start, ranges[i].size,
                      ranges[j].size) +=
          smallmatrix.block(inv_ranges[i].start, inv_ranges[j].start,
                            inv_ranges[i].size, inv_ranges[j].size);
    }
  }
  return;
}

Eigen::MatrixXd GridBox::ReadFromBigMatrix(
    const Eigen::MatrixXd& bigmatrix) const {
  Eigen::MatrixXd matrix = Eigen::MatrixXd(matrix_size, matrix_size);
  for (unsigned i = 0; i < ranges.size(); i++) {
    for (unsigned j = 0; j < ranges.size(); j++) {
      matrix.block(inv_ranges[i].start, inv_ranges[j].start, inv_ranges[i].size,
                   inv_ranges[j].size) =
          bigmatrix.block(ranges[i].start, ranges[j].start, ranges[i].size,
                          ranges[j].size);
    }
  }
  return matrix;
}

void GridBox::PrepareForIntegration() {
  unsigned index = 0;
  aoranges = std::vector<GridboxRange>(0);
  ranges = std::vector<GridboxRange>(0);
  inv_ranges = std::vector<GridboxRange>(0);
  std::vector<unsigned> start;
  std::vector<unsigned> end;

  for (const auto shell : significant_shells) {
    GridboxRange temp;
    temp.size = shell->getNumFunc();
    temp.start = index;
    aoranges.push_back(temp);
    index += shell->getNumFunc();
    start.push_back(shell->getStartIndex());
    end.push_back(shell->getStartIndex() + shell->getNumFunc());
  }
  std::vector<unsigned> startindex;
  std::vector<unsigned> endindex;

  if (start.size() > 1) {
    startindex.push_back(start[0]);

    for (unsigned i = 0; i < start.size() - 1; ++i) {

      if (end[i] != start[i + 1]) {
        startindex.push_back(start[i + 1]);
        endindex.push_back(end[i]);
      }
    }
    endindex.push_back(end[end.size() - 1]);
  } else {
    startindex = start;
    endindex = end;
  }
  unsigned shellstart = 0;
  for (unsigned i = 0; i < startindex.size(); ++i) {
    unsigned size = endindex[i] - startindex[i];
    GridboxRange temp;
    temp.size = size;
    temp.start = startindex[i];
    ranges.push_back(temp);
    GridboxRange temp2;
    temp2.size = size;
    temp2.start = shellstart;
    inv_ranges.push_back(temp2);
    shellstart += size;
  }
  return;
}

}  // namespace xtp
}  // namespace votca