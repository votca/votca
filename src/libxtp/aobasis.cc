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

// Standard includes
#include <vector>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/qmpackage.h"

namespace votca {
namespace xtp {

Index AOBasis::getMaxL() const {
  Index n = 0;
  for (const auto& shell : _aoshells) {
    n = std::max(static_cast<Index>(shell.getL()), n);
  }
  return n;
}

Index AOBasis::getMaxNprim() const {
  Index n = 0;
  for (const auto& shell : _aoshells) {
    n = std::max(shell.getSize(), n);
  }
  return n;
}

std::vector<Index> AOBasis::getMapToBasisFunctions() const {
  std::vector<Index> result;
  result.reserve(_aoshells.size());

  Index n = 0;
  for (const auto& shell : _aoshells) {
    result.push_back(n);
    n += shell.getNumFunc();
  }
  return result;
}

AOShell& AOBasis::addShell(const Shell& shell, const QMAtom& atom,
                           Index startIndex) {
  _aoshells.push_back(AOShell(shell, atom, startIndex));
  return _aoshells.back();
}

const std::vector<const AOShell*> AOBasis::getShellsofAtom(Index AtomId) const {
  std::vector<const AOShell*> result;
  for (const auto& aoshell : _aoshells) {
    if (aoshell.getAtomIndex() == AtomId) {
      result.push_back(&aoshell);
    }
  }
  return result;
}

void AOBasis::Fill(const BasisSet& bs, const QMMolecule& atoms) {
  _name = bs.Name();
  _AOBasisSize = 0;
  _aoshells.clear();
  _FuncperAtom.clear();
  // loop over atoms
  for (const QMAtom& atom : atoms) {
    Index atomfunc = 0;
    const std::string& name = atom.getElement();
    const Element& element = bs.getElement(name);
    for (const Shell& shell : element) {
      Index numfuncshell = NumFuncShell(shell.getL());
      AOShell& aoshell = addShell(shell, atom, _AOBasisSize);
      _AOBasisSize += numfuncshell;
      atomfunc += numfuncshell;
      for (const GaussianPrimitive& gaussian : shell) {
        aoshell.addGaussian(gaussian);
      }
      aoshell.CalcMinDecay();
      aoshell.normalizeContraction();
    }
    _FuncperAtom.push_back(atomfunc);
  }
  GenerateLibintBasis();
  return;
}

std::vector<libint2::Shell> AOBasis::GenerateLibintBasis() const {
  std::vector<libint2::Shell> libintshells;
  libintshells.reserve(_aoshells.size());
  for (const auto& shell : _aoshells) {
    libintshells.push_back(shell.LibintShell());
  }
  return libintshells;
}

}  // namespace xtp
}  // namespace votca
