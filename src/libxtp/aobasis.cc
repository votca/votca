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
#include "votca/xtp/aoshell.h"
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
  clear();
  _name = bs.Name();
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
  }
  FillFuncperAtom();
  GenerateLibintBasis();
  return;
}

void AOBasis::FillFuncperAtom() {
  _FuncperAtom.clear();
  Index currentIndex = -1;
  for (const auto& shell : _aoshells) {
    if (shell.getAtomIndex() == currentIndex) {
      _FuncperAtom[shell.getAtomIndex()] += shell.getNumFunc();
    } else {
      currentIndex = shell.getAtomIndex();
      _FuncperAtom.push_back(shell.getNumFunc());
    }
  }
}

void AOBasis::GenerateLibintBasis() {
  _libintshells.clear();
  _libintshells.reserve(_aoshells.size());

  for (const auto& shell : _aoshells) {
    libint2::svector<libint2::Shell::real_t> decays;
    libint2::svector<libint2::Shell::Contraction> contractions;
    const Eigen::Vector3d& pos = shell.getPos();
    libint2::Shell::Contraction contr;
    contr.l = static_cast<int>(shell.getL());
    contr.pure = true;
    for (const auto& primitive : shell) {
      decays.push_back(primitive.getDecay());
      contr.coeff.push_back(primitive.getContraction());
    }
    contractions.push_back(contr);
    std::array<libint2::Shell::real_t, 3> libintpos = {pos[0], pos[1], pos[2]};
    libint2::Shell libintshell(decays, contractions, libintpos);
    _libintshells.push_back(libintshell);
  }
}

void AOBasis::clear() {
  _name = "";
  _aoshells.clear();
  _libintshells.clear();
  _FuncperAtom.clear();
  _AOBasisSize = 0;
}

void AOBasis::WriteToCpt(CheckpointWriter& w) const {
  w(_name, "name");

  Index numofprimitives = 0;
  for (const auto& shell : _aoshells) {
    numofprimitives += shell.getSize();
  }

  const AOGaussianPrimitive& dummy = *(getShell(0).begin());
  CptTable table = w.openTable("Contractions", dummy, numofprimitives);

  std::vector<AOGaussianPrimitive::data> dataVec(numofprimitives);
  Index i = 0;
  for (const auto& shell : _aoshells) {
    for (const auto& gaussian : shell) {
      gaussian.WriteData(dataVec[i]);
      i++;
    }
  }

  table.write(dataVec);
}

void AOBasis::ReadFromCpt(CheckpointReader& r) {
  clear();
  r(_name, "name");

  FillFuncperAtom();
  GenerateLibintBasis();
}

}  // namespace xtp
}  // namespace votca
