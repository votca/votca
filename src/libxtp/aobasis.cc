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

void AOBasis::GenerateLibintBasis() {
  _libintshells.resize(0);
  _libintshells.reserve(_aoshells.size());

  for (const auto& shell : _aoshells) {
    std::vector<libint2::real_t> decays;
    std::vector<libint2::Shell::Contraction> contractions;
    const Eigen::Vector3d& pos = shell.getPos();
    for (const auto& primitive : shell) {
      decays.push_back(primitive.getDecay());
      libint2::Shell::Contraction contr;
      contr.l = static_cast<int>(shell.getL());
      contr.pure = true;
      contr.coeff.push_back(primitive.getContraction());
      contractions.push_back(contr);
    }

    std::array<libint2::real_t, 3> libintpos = {pos[0], pos[1], pos[2]};
    libint2::Shell libintshell(decays, contractions, libintpos);
    _libintshells.push_back(libintshell);
  }
}

}  // namespace xtp
}  // namespace votca
