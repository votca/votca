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

// Local VOTCA includes
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/basisset.h"

namespace votca {
namespace xtp {

std::vector<Transposition> OrbReorder::computeTranspositions(
    std::vector<Index> vStart, std::vector<Index> vTarget) const {
  if (vStart.size() != vTarget.size()) {
    throw std::runtime_error(
        "Can't compute transpositions, reorder vectors have different "
        "size!\n");
  }
  std::vector<Transposition> transpositions;
  for (Index i = 0; i < static_cast<Index>(vStart.size()); i++) {
    std::vector<Index>::iterator it;
    it = std::find(vStart.begin(), vStart.end(), vTarget[i]);
    Index pos = std::distance(vStart.begin(), it);
    std::swap(vStart[i], vStart[pos]);
    if (i != pos) {
      transpositions.push_back(Transposition{i, pos});
    }
  }
  return transpositions;
}

std::vector<Index> OrbReorder::copySegment(const std::array<Index, 49>& input,
                                           Index start, Index size) const {
  return std::vector<Index>{input.begin() + start,
                            input.begin() + start + size};
}

OrbReorder::OrbReorder(std::array<Index, 49> reorder,
                       std::array<Index, 49> multipliers, bool reverse)
    : multipliers_(multipliers), reorder_(reorder), reverse_(reverse) {

  // Compute transpositions for every individual shell
  Index currentFunction = 0;
  for (int l = 0; l < 7; l++) {
    Index nrOfFunctions = NumFuncShell(static_cast<L>(l));
    if (!reverse_) {  // ordering from external to votca order
      transpositions_[l] = computeTranspositions(
          copySegment(reorder_, currentFunction, nrOfFunctions),
          copySegment(votcaOrder_, currentFunction, nrOfFunctions));
    } else {  // votca order to external order
      transpositions_[l] = computeTranspositions(
          copySegment(votcaOrder_, currentFunction, nrOfFunctions),
          copySegment(reorder_, currentFunction, nrOfFunctions));
    }
    currentFunction += nrOfFunctions;
  }
}

void OrbReorder::reorderOrbitals(Eigen::MatrixXd& moCoefficients,
                                 const AOBasis& basis) {

  for (const AOShell& shell : basis) {
    Index currentFunction = shell.getStartIndex();

    if (reverse_) {  // multiply first before reversing ordering
      // Get multiplier vector for shell
      std::vector<Index> shellmultiplier =
          copySegment(multipliers_, shell.getOffset(), shell.getNumFunc());

      // multiply shell
      for (Index i = 0; i < shell.getNumFunc(); i++) {
        moCoefficients.row(currentFunction + i) *= double(shellmultiplier[i]);
      }
    }

    // reorder shell
    Index l = static_cast<Index>(shell.getL());
    for (const Transposition& transposition : transpositions_[l]) {
      moCoefficients.row(currentFunction + transposition.first)
          .swap(moCoefficients.row(currentFunction + transposition.second));
    }

    if (!reverse_) {  // multiply after reordering
      // Get multiplier vector for shell
      std::vector<Index> shellmultiplier =
          copySegment(multipliers_, shell.getOffset(), shell.getNumFunc());

      // multiply shell
      for (Index i = 0; i < shell.getNumFunc(); i++) {
        moCoefficients.row(currentFunction + i) *= double(shellmultiplier[i]);
      }
    }
  }
}
void OrbReorder::reorderOperator(Eigen::MatrixXd& Matrixoperator,
                                 const AOBasis& basis) {
  // reorder rows first
  reorderOrbitals(Matrixoperator, basis);
  Matrixoperator.transposeInPlace();
  reorderOrbitals(Matrixoperator, basis);
  Matrixoperator.transposeInPlace();
}

}  // namespace xtp
}  // namespace votca