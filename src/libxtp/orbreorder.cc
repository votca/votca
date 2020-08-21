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

namespace votca {
namespace xtp {

std::vector<std::array<Index, 2>> OrbReorder::computeTranspositions(
    std::vector<Index> vStart, std::vector<Index> vTarget) {
  if (vStart.size() != vTarget.size()) {
    throw "Can't compute transpositions, reorder vectors have different "
          "size!\n";
  }
  std::vector<std::array<Index, 2>> transpositions;
  for (Index i = 0; i < static_cast<Index>(vStart.size()); i++) {
    std::vector<Index>::iterator it;
    it = std::find(vStart.begin(), vStart.end(), vTarget[i]);
    Index pos = std::distance(vStart.begin(), it);
    std::swap(vStart[i], vStart[pos]);
    if (i != pos) {
      transpositions.push_back(std::array<Index, 2>{i, pos});
    }
  }
  return transpositions;
}

OrbReorder::OrbReorder(std::array<Index, 25> reorder,
                       std::array<Index, 25> multipliers, bool reverse)
    : _multipliers(multipliers), _reorder(reorder) {

  // Compute transpositions for every shell
  int currentFunction = 0;
  for (int l = 0; l < 5; l++) {
    int nrOfFunctions = 2 * l + 1;
    if (!reverse) {
      _transpositions[l] = computeTranspositions(
          std::vector<Index>{
              _reorder.begin() + currentFunction,
              _reorder.begin() + currentFunction + nrOfFunctions},
          std::vector<Index>{
              _votcaOrder.begin() + currentFunction,
              _votcaOrder.begin() + currentFunction + nrOfFunctions});
    } else {
      _transpositions[l] = computeTranspositions(
          std::vector<Index>{
              _votcaOrder.begin() + currentFunction,
              _votcaOrder.begin() + currentFunction + nrOfFunctions},
          std::vector<Index>{
              _reorder.begin() + currentFunction,
              _reorder.begin() + currentFunction + nrOfFunctions});
    }
    for (auto& tr : _transpositions[l]) {
      std::cout << tr[0] << " " << tr[1] << std::endl;
    }
    currentFunction += nrOfFunctions;
  }
}

void OrbReorder::reorderOrbitals(Eigen::MatrixXd& moCoefficients,
                                 const AOBasis& basis) {
  std::vector<Index> multiplier;
  multiplier.reserve(basis.AOBasisSize());

  // reorder and get multiplier vector
  Index currentFunction = 0;
  for (const AOShell& shell : basis) {
    // Make multiplier vector
    std::vector<Index> shellmultiplier{
        _multipliers.begin() + shell.getOffset(),
        _multipliers.begin() + shell.getOffset() + shell.getNumFunc()};
    multiplier.insert(multiplier.end(), shellmultiplier.begin(),
                      shellmultiplier.end());
    // reorder
    Index l = static_cast<Index>(shell.getL());
    for (auto& transposition : _transpositions[l]) {
      moCoefficients.row(currentFunction + transposition[0])
          .swap(moCoefficients.row(currentFunction + transposition[1]));
    }
    currentFunction += (2 * l + 1);
  }

  // Multiply by multiplier
  for (Index i = 0; i < moCoefficients.cols(); i++) {
    moCoefficients.row(i) = multiplier[i] * moCoefficients.row(i);
  }
}  // namespace xtp

}  // namespace xtp
}  // namespace votca