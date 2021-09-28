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

#pragma once
#ifndef VOTCA_XTP_AOBASIS_H
#define VOTCA_XTP_AOBASIS_H

// Local VOTCA includes
#include "aoshell.h"

#include "eigen.h"
#include <libint2/shell.h>

namespace votca {
namespace xtp {
class QMMolecule;
class BasisSet;
class CheckpointWriter;
class CheckpointReader;

/**
 * \brief Container to hold Basisfunctions for all atoms
 *
 * It is constructed from a QMMolecule and a BasisSet.
 */
class AOBasis {
 public:
  void Fill(const BasisSet& bs, const QMMolecule& atoms);

  Index AOBasisSize() const { return AOBasisSize_; }

  using AOShellIterator = std::vector<AOShell>::const_iterator;
  AOShellIterator begin() const { return aoshells_.begin(); }
  AOShellIterator end() const { return aoshells_.end(); }

  const AOShell& getShell(Index idx) const { return aoshells_[idx]; }

  const std::vector<const AOShell*> getShellsofAtom(Index AtomId) const;

  Index getNumofShells() const { return Index(aoshells_.size()); }

  Index getNumberOfPrimitives() const {
    Index totalPrimitives = 0;
    for (const AOShell& shell : aoshells_) {
      totalPrimitives += shell.getSize();
    }
    return totalPrimitives;
  }

  Index getMaxNprim() const;

  Index getMaxL() const;

  std::vector<Index> getMapToBasisFunctions() const;

  const std::vector<Index>& getFuncPerAtom() const { return FuncperAtom_; }

  std::vector<libint2::Shell> GenerateLibintBasis() const;

  std::vector<std::vector<Index>> ComputeShellPairs(
      double threshold = 1e-20) const;

  AOShell& addShell(const Shell& shell, const QMAtom& atom, Index startIndex);

  const std::string& Name() const { return name_; }

  void UpdateShellPositions(const QMMolecule& mol);

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r);

  void add(const AOBasis& other);

  friend std::ostream& operator<<(std::ostream& out, const AOBasis& aobasis);

 private:
  void FillFuncperAtom();

  void clear();
  std::string name_ = "";

  std::vector<AOShell> aoshells_;

  std::vector<Index> FuncperAtom_;

  Index AOBasisSize_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOBASIS_H
