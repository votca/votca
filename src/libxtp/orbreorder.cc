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

std::vector<Index> OrbReorder::copySegment(const std::array<Index, 25>& input,
                                           Index start, Index size) const {
  return std::vector<Index>{input.begin() + start,
                            input.begin() + start + size};
}

OrbReorder::OrbReorder(std::array<Index, 25> reorder,
                       std::array<Index, 25> multipliers, bool reverse)
    : _multipliers(multipliers), _reorder(reorder) {

  // Compute transpositions for every individual shell
  Index currentFunction = 0;
  for (int l = 0; l < 5; l++) {
    Index nrOfFunctions = NumFuncShell(static_cast<L>(l));
    if (!reverse) {
      _transpositions[l] = computeTranspositions(
          copySegment(_reorder, currentFunction, nrOfFunctions),
          copySegment(_votcaOrder, currentFunction, nrOfFunctions));
    } else {
      _transpositions[l] = computeTranspositions(
          copySegment(_votcaOrder, currentFunction, nrOfFunctions),
          copySegment(_reorder, currentFunction, nrOfFunctions));
    }
    currentFunction += nrOfFunctions;
  }
}  // namespace xtp

void OrbReorder::reorderOrbitals(Eigen::MatrixXd& moCoefficients,
                                 const AOBasis& basis) {
  std::vector<Index> multiplier;
  multiplier.reserve(basis.AOBasisSize());

  Index currentFunction = 0;
  for (const AOShell& shell : basis) {

    // reorder shell
    Index l = static_cast<Index>(shell.getL());
    for (const Transposition& transposition : _transpositions[l]) {
      moCoefficients.row(currentFunction + transposition.first)
          .swap(moCoefficients.row(currentFunction + transposition.second));
    }

    // Get multiplier vector for shell
    std::vector<Index> shellmultiplier =
        copySegment(_multipliers, shell.getOffset(), shell.getNumFunc());

    // multiply shell
    for (Index i = 0; i < shell.getNumFunc(); i++) {
      moCoefficients.row(currentFunction + i) *= double(shellmultiplier[i]);
    }
    currentFunction += shell.getNumFunc();
  }
}
void OrbReorder::reorderRowsAndCols(Eigen::MatrixXd& moCoefficients,
                                    const AOBasis& basis) {
  // reorder rows first
  reorderOrbitals(moCoefficients, basis);

  // next the cols
  Index currentFunction = 0;
  for (const AOShell& shell : basis) {

    // reorder shell
    Index l = static_cast<Index>(shell.getL());
    for (const Transposition& transposition : _transpositions[l]) {
      moCoefficients.col(currentFunction + transposition.first)
          .swap(moCoefficients.col(currentFunction + transposition.second));
    }
    currentFunction += shell.getNumFunc();
  }
}

}  // namespace xtp
}  // namespace votca