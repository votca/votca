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

#include <votca/xtp/gridbox.h>

namespace votca {
namespace xtp {

void GridBox::FindSignificantShells(const AOBasis& basis) {

  for (const AOShell& store : basis) {
    const double decay = store.getMinDecay();
    const Eigen::Vector3d& shellpos = store.getPos();
    for (const auto& point : grid_pos) {
      Eigen::Vector3d dist = shellpos - point;
      double distsq = dist.squaredNorm();
      // if contribution is smaller than -ln(1e-10), add shell to list
      if ((decay * distsq) < 20.7) {
        addShell(&store);
        break;
      }
    }
  }
}

Eigen::VectorXd GridBox::CalcAOValue_and_Grad(
    Eigen::MatrixX3d& ao_grad, const Eigen::Vector3d& point) const {
  Eigen::VectorXd ao = Eigen::VectorXd::Zero(Matrixsize());
  for (Index j = 0; j < Shellsize(); ++j) {
    Eigen::Block<Eigen::MatrixX3d> grad_block =
        ao_grad.block(aoranges[j].start, 0, aoranges[j].size, 3);
    Eigen::VectorBlock<Eigen::VectorXd> ao_block =
        ao.segment(aoranges[j].start, aoranges[j].size);
    significant_shells[j]->EvalAOspace(ao_block, grad_block, point);
  }
  return ao;
}

Eigen::VectorXd GridBox::CalcAOValues(const Eigen::Vector3d& pos) const {
  Eigen::VectorXd ao = Eigen::VectorXd::Zero(Matrixsize());
  for (Index j = 0; j < Shellsize(); ++j) {
    Eigen::VectorBlock<Eigen::VectorXd> ao_block =
        ao.segment(aoranges[j].start, aoranges[j].size);
    significant_shells[j]->EvalAOspace(ao_block, pos);
  }
  return ao;
}

void GridBox::AddtoBigMatrix(Eigen::MatrixXd& bigmatrix,
                             const Eigen::MatrixXd& smallmatrix) const {
  for (Index i = 0; i < Index(ranges.size()); i++) {
    for (Index j = 0; j < Index(ranges.size()); j++) {
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
  for (Index i = 0; i < Index(ranges.size()); i++) {
    for (Index j = 0; j < Index(ranges.size()); j++) {
      matrix.block(inv_ranges[i].start, inv_ranges[j].start, inv_ranges[i].size,
                   inv_ranges[j].size) =
          bigmatrix.block(ranges[i].start, ranges[j].start, ranges[i].size,
                          ranges[j].size);
    }
  }
  return matrix;
}

void GridBox::PrepareForIntegration() {
  Index index = 0;
  aoranges = std::vector<GridboxRange>(0);
  ranges = std::vector<GridboxRange>(0);
  inv_ranges = std::vector<GridboxRange>(0);
  std::vector<Index> start;
  std::vector<Index> end;

  for (const auto shell : significant_shells) {
    GridboxRange temp;
    temp.size = shell->getNumFunc();
    temp.start = index;
    aoranges.push_back(temp);
    index += shell->getNumFunc();
    start.push_back(shell->getStartIndex());
    end.push_back(shell->getStartIndex() + shell->getNumFunc());
  }
  std::vector<Index> startindex;
  std::vector<Index> endindex;

  if (start.size() > 1) {
    startindex.push_back(start[0]);

    for (Index i = 0; i < Index(start.size()) - 1; ++i) {

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
  Index shellstart = 0;
  for (Index i = 0; i < Index(startindex.size()); ++i) {
    Index size = endindex[i] - startindex[i];
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