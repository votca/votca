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

#pragma once
#ifndef VOTCA_XTP_ECPAOBASIS_H
#define VOTCA_XTP_ECPAOBASIS_H
// Local VOTCA includes
#include "checkpoint.h"
#include "qmatom.h"
#include "ecpbasisset.h"
#include <libecpint/ecp.hpp>

namespace votca {
namespace xtp {
class QMMolecule;
class ECPBasisSet;


 std::ostream& operator<<(std::ostream& out, const  libecpint::ECP& ecp);


/**
 * \brief Container to hold ECPs for all atoms
 *
 * It is constructed from a vector of QMAtoms and an ECPBasisSet.
 */
class ECPAOBasis {
 public:
  // returns element names for which no ecp was found
  std::vector<std::string> Fill(const ECPBasisSet& bs, QMMolecule& atoms);

  Index ECPAOBasisSize() const { return _AOBasisSize; }

  using constECPAOShellIterator = std::vector<libecpint::ECP>::const_iterator;
  constECPAOShellIterator begin() const { return _aoshells.begin(); }
  constECPAOShellIterator end() const { return _aoshells.end(); }

  using ECPAOShellIterator = std::vector<libecpint::ECP>::iterator;
  ECPAOShellIterator begin()  { return _aoshells.begin(); }
  ECPAOShellIterator end() { return _aoshells.end(); }

  Index getMaxL() const;
  void AddECPChargeToMolecule(QMMolecule& mol) const;

  const std::string& Name() const { return _name; }

  void UpdateShellPositions(const QMMolecule& mol);

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r);

  void add(const ECPAOBasis& other);

  friend std::ostream& operator<<(std::ostream& out, const ECPAOBasis& ecp);

 private:
  void clear();

  libecpint::ECP& addShell(const ECPShell& shell, const QMAtom& atom);

  std::vector<Index> _ncore_perAtom;

  std::vector<libecpint::ECP> _aoshells;

  std::string _name = "";
  Index _AOBasisSize;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ECPAOBASIS_H
