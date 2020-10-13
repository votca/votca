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
#include "votca/xtp/aobasis.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/make_libint_work.h"
#include "votca/xtp/qmmolecule.h"
#include <libint2.hpp>

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

void AOBasis::add(const AOBasis& other) {
  Index atomindex_offset = Index(_FuncperAtom.size());
  for (AOShell shell : other) {
    shell._atomindex += atomindex_offset;
    shell._startIndex = _AOBasisSize;
    _AOBasisSize += shell.getNumFunc();
    _aoshells.push_back(shell);
  }

  FillFuncperAtom();
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
std::vector<libint2::Shell> AOBasis::GenerateLibintBasis() const {
  std::vector<libint2::Shell> libintshells;
  libintshells.reserve(_aoshells.size());
  for (const auto& shell : _aoshells) {
    libintshells.push_back(shell.LibintShell());
  }
  return libintshells;
}

std::vector<std::vector<Index>> AOBasis::ComputeShellPairs(
    double threshold) const {

  Index nthreads = OPENMP::getMaxThreads();

  std::vector<libint2::Shell> shells = GenerateLibintBasis();

  // construct the 2-electron repulsion integrals engine
  std::vector<libint2::Engine> engines;
  engines.reserve(nthreads);
  engines.emplace_back(libint2::Operator::overlap, getMaxNprim(), getMaxL(), 0);
  for (Index i = 1; i != nthreads; ++i) {
    engines.push_back(engines[0]);
  }

  std::vector<std::vector<Index>> pairs(shells.size());

#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0; s1 < Index(shells.size()); ++s1) {
    Index thread_id = OPENMP::getThreadId();

    libint2::Engine& engine = engines[thread_id];
    const libint2::Engine::target_ptr_vec& buf = engine.results();
    Index n1 = shells[s1].size();

    for (Index s2 = s1; s2 < Index(shells.size()); ++s2) {
      bool on_same_center = (shells[s1].O == shells[s2].O);
      bool significant = on_same_center;
      if (!on_same_center) {
        Index n2 = shells[s2].size();
        engine.compute(shells[s1], shells[s2]);
        Eigen::Map<const Eigen::MatrixXd> buf_mat(buf[0], n1, n2);
        significant = (buf_mat.norm() >= threshold);
      }
      if (significant) {
        pairs[s1].push_back(s2);
      }
    }
    std::vector<Index>& list = pairs[s1];
    std::sort(list.begin(), list.end());
  }
  return pairs;
}

void AOBasis::UpdateShellPositions(const QMMolecule& mol) {
  for (AOShell& shell : _aoshells) {
    shell._pos = mol[shell.getAtomIndex()].getPos();
  }
}

void AOBasis::clear() {
  _name = "";
  _aoshells.clear();
  _FuncperAtom.clear();
  _AOBasisSize = 0;
}

void AOBasis::WriteToCpt(CheckpointWriter& w) const {
  w(_name, "name");
  w(_AOBasisSize, "basissize");
  Index numofprimitives = 0;
  for (const auto& shell : _aoshells) {
    numofprimitives += shell.getSize();
  }

  // this is all to make dummy AOGaussian
  Shell s(L::S, 0);
  GaussianPrimitive d(0.1, 0.1);
  QMAtom dummy(0, "H", Eigen::Vector3d::Zero());
  AOShell s1(s, dummy, 0);
  s1.addGaussian(d);
  const AOGaussianPrimitive& dummy2 = *s1.begin();

  CptTable table = w.openTable("Contractions", dummy2, numofprimitives);

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
  r(_AOBasisSize, "basissize");
  if (_AOBasisSize > 0) {
    // this is all to make dummy AOGaussian
    Shell s(L::S, 0);
    GaussianPrimitive d(0.1, 0.1);
    QMAtom dummy(0, "H", Eigen::Vector3d::Zero());
    AOShell s1(s, dummy, 0);
    s1.addGaussian(d);
    const AOGaussianPrimitive& dummy2 = *s1.begin();

    CptTable table = r.openTable("Contractions", dummy2);
    std::vector<AOGaussianPrimitive::data> dataVec(table.numRows());
    table.read(dataVec);
    Index laststartindex = -1;
    for (std::size_t i = 0; i < table.numRows(); ++i) {
      if (dataVec[i].startindex != laststartindex) {
        _aoshells.push_back(AOShell(dataVec[i]));
        laststartindex = dataVec[i].startindex;
      } else {
        _aoshells.back()._gaussians.push_back(
            AOGaussianPrimitive(dataVec[i], _aoshells.back()));
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
