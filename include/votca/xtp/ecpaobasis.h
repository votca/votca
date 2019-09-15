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
#ifndef VOTCA_XTP_ECPAOBasis_H
#define VOTCA_XTP_ECPAOBasis_H

#include <votca/xtp/ecpaoshell.h>
#include <votca/xtp/eigen.h>

namespace votca {
namespace xtp {
class QMMolecule;
class ECPBasisSet;

/**
 * \brief Container to hold ECPs for all atoms
 *
 * It is constructed from a vector of QMAtoms and an ECPBasisSet.
 */
class ECPAOBasis {
 public:
  // returns element names for which no ecp was found
  std::vector<std::string> Fill(const ECPBasisSet& bs, QMMolecule& atoms);

  int ECPAOBasisSize() const { return _AOBasisSize; }

  typedef std::vector<ECPAOShell>::const_iterator ECPAOShellIterator;
  ECPAOShellIterator begin() const { return _aoshells.begin(); }
  ECPAOShellIterator end() const { return _aoshells.end(); }

  const ECPAOShell& getShell(int idx) const { return _aoshells[idx]; }

  const ECPAOShell& back() const { return _aoshells.back(); }

  const std::vector<std::vector<const ECPAOShell*> >& ShellsPerAtom() const {
    return _shells_perAtom;
  }

 private:
  ECPAOShell& addShell(const ECPShell& shell, const QMAtom& atom,
                       int startIndex, int Lmax);

  std::vector<ECPAOShell> _aoshells;

  std::vector<std::vector<const ECPAOShell*> > _shells_perAtom;
  int _AOBasisSize;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ECPAOBasis_H
