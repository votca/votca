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

#pragma once
#ifndef VOTCA_XTP_AOBASIS_H
#define VOTCA_XTP_AOBASIS_H

#include <votca/xtp/aoshell.h>
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {
class QMMolecule;
class BasisSet;

/**
 * \brief Container to hold Basisfunctions for all atoms
 *
 * It is constructed from a QMMolecule and a BasisSet.
 */
class AOBasis {
 public:
  void ReorderMOs(Eigen::MatrixXd& v, const std::string& start,
                  const std::string& target) const;

  void Fill(const BasisSet& bs, const QMMolecule& atoms);

  Index AOBasisSize() const { return _AOBasisSize; }

  using AOShellIterator = std::vector<AOShell>::const_iterator;
  AOShellIterator begin() const { return _aoshells.begin(); }
  AOShellIterator end() const { return _aoshells.end(); }

  const AOShell& getShell(Index idx) const { return _aoshells[idx]; }

  const std::vector<const AOShell*> getShellsofAtom(Index AtomId) const;

  Index getNumofShells() const { return Index(_aoshells.size()); }

  const std::vector<Index>& getFuncPerAtom() const { return _FuncperAtom; }

 private:
  AOShell& addShell(const Shell& shell, const QMAtom& atom, Index startIndex);

  void MultiplyMOs(Eigen::MatrixXd& v,
                   const std::vector<Index>& multiplier) const;

  std::vector<Index> invertOrder(const std::vector<Index>& order) const;

  std::vector<Index> getReorderVector(const std::string& start,
                                      const std::string& target) const;

  void addReorderShell(const std::string& start, const std::string& target,
                       const std::string& shell,
                       std::vector<Index>& neworder) const;

  std::vector<Index> getMultiplierVector(const std::string& start,
                                         const std::string& target) const;

  void addMultiplierShell(const std::string& start, const std::string& target,
                          const std::string& shell,
                          std::vector<Index>& multiplier) const;

  std::vector<AOShell> _aoshells;

  std::vector<Index> _FuncperAtom;

  Index _AOBasisSize;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOBASIS_H