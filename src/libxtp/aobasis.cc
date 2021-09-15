/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
#include "votca/xtp/aobasis.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/qmmolecule.h"

namespace votca {
namespace xtp {

Index AOBasis::getMaxL() const {
  Index n = 0;
  for (const auto& shell : aoshells_) {
    n = std::max(static_cast<Index>(shell.getL()), n);
  }
  return n;
}

Index AOBasis::getMaxNprim() const {
  Index n = 0;
  for (const auto& shell : aoshells_) {
    n = std::max(shell.getSize(), n);
  }
  return n;
}

std::vector<Index> AOBasis::getMapToBasisFunctions() const {
  std::vector<Index> result;
  result.reserve(aoshells_.size());

  Index n = 0;
  for (const auto& shell : aoshells_) {
    result.push_back(n);
    n += shell.getNumFunc();
  }
  return result;
}

AOShell& AOBasis::addShell(const Shell& shell, const QMAtom& atom,
                           Index startIndex) {
  aoshells_.push_back(AOShell(shell, atom, startIndex));
  return aoshells_.back();
}

const std::vector<const AOShell*> AOBasis::getShellsofAtom(Index AtomId) const {
  std::vector<const AOShell*> result;
  for (const auto& aoshell : aoshells_) {
    if (aoshell.getAtomIndex() == AtomId) {
      result.push_back(&aoshell);
    }
  }
  return result;
}

void AOBasis::add(const AOBasis& other) {
  Index atomindex_offset = Index(FuncperAtom_.size());
  for (AOShell shell : other) {
    shell.atomindex_ += atomindex_offset;
    shell.startIndex_ = AOBasisSize_;
    AOBasisSize_ += shell.getNumFunc();
    aoshells_.push_back(shell);
  }

  FillFuncperAtom();
}

void AOBasis::Fill(const BasisSet& bs, const QMMolecule& atoms) {
  clear();
  name_ = bs.Name();
  // loop over atoms
  for (const QMAtom& atom : atoms) {
    Index atomfunc = 0;
    const std::string& name = atom.getElement();
    const Element& element = bs.getElement(name);
    for (const Shell& shell : element) {
      Index numfuncshell = NumFuncShell(shell.getL());
      AOShell& aoshell = addShell(shell, atom, AOBasisSize_);
      AOBasisSize_ += numfuncshell;
      atomfunc += numfuncshell;
      for (const GaussianPrimitive& gaussian : shell) {
        aoshell.addGaussian(gaussian);
      }
      aoshell.CalcMinDecay();
      aoshell.normalizeContraction();
    }
  }
  FillFuncperAtom();
  return;
}

void AOBasis::FillFuncperAtom() {
  FuncperAtom_.clear();
  Index currentIndex = -1;
  for (const auto& shell : aoshells_) {
    if (shell.getAtomIndex() == currentIndex) {
      FuncperAtom_[shell.getAtomIndex()] += shell.getNumFunc();
    } else {
      currentIndex = shell.getAtomIndex();
      FuncperAtom_.push_back(shell.getNumFunc());
    }
  }
}
std::vector<libint2::Shell> AOBasis::GenerateLibintBasis() const {
  std::vector<libint2::Shell> libintshells;
  libintshells.reserve(aoshells_.size());
  for (const auto& shell : aoshells_) {
    libintshells.push_back(shell.LibintShell());
  }
  return libintshells;
}

void AOBasis::UpdateShellPositions(const QMMolecule& mol) {
  for (AOShell& shell : aoshells_) {
    shell.pos_ = mol[shell.getAtomIndex()].getPos();
  }
}

void AOBasis::clear() {
  name_ = "";
  aoshells_.clear();
  FuncperAtom_.clear();
  AOBasisSize_ = 0;
}

void AOBasis::WriteToCpt(CheckpointWriter& w) const {
  w(name_, "name");
  w(AOBasisSize_, "basissize");
  Index numofprimitives = 0;
  for (const auto& shell : aoshells_) {
    numofprimitives += shell.getSize();
  }

  CptTable table =
      w.openTable<AOGaussianPrimitive>("Contractions", numofprimitives);

  std::vector<AOGaussianPrimitive::data> dataVec(numofprimitives);
  Index i = 0;
  for (const auto& shell : aoshells_) {
    for (const auto& gaussian : shell) {
      gaussian.WriteData(dataVec[i], shell);
      i++;
    }
  }

  table.write(dataVec);
}

void AOBasis::ReadFromCpt(CheckpointReader& r) {
  clear();
  r(name_, "name");
  r(AOBasisSize_, "basissize");
  if (AOBasisSize_ > 0) {

    CptTable table = r.openTable<AOGaussianPrimitive>("Contractions");
    std::vector<AOGaussianPrimitive::data> dataVec(table.numRows());
    table.read(dataVec);
    Index laststartindex = -1;
    for (std::size_t i = 0; i < table.numRows(); ++i) {
      if (dataVec[i].startindex != laststartindex) {
        aoshells_.push_back(AOShell(dataVec[i]));
        laststartindex = dataVec[i].startindex;
      } else {
        aoshells_.back().gaussians_.push_back(AOGaussianPrimitive(dataVec[i]));
      }
    }

    FillFuncperAtom();
  }
}

std::ostream& operator<<(std::ostream& out, const AOBasis& aobasis) {
  out << "Name:" << aobasis.Name() << "\n";
  out << " Functions:" << aobasis.AOBasisSize()
      << " Shells:" << aobasis.getNumofShells() << "\n";
  for (const auto& shell : aobasis) {
    out << shell;
  }
  return out;
}

}  // namespace xtp
}  // namespace votca
