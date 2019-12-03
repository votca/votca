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

#include <votca/xtp/aomatrix.h>
#include <votca/xtp/couplingbase.h>

namespace votca {
namespace xtp {

Eigen::MatrixXd CouplingBase::CalculateOverlapMatrix(
    const Orbitals& orbitalsAB) const {
  AOBasis dftbasis = orbitalsAB.SetupDftBasis();
  AOOverlap dftAOoverlap;
  dftAOoverlap.Fill(dftbasis);
  return dftAOoverlap.Matrix();
}

void CouplingBase::CheckAtomCoordinates(const Orbitals& orbitalsA,
                                        const Orbitals& orbitalsB,
                                        const Orbitals& orbitalsAB) const {
  const QMMolecule& atomsA = orbitalsA.QMAtoms();
  const QMMolecule& atomsB = orbitalsB.QMAtoms();
  const QMMolecule& atomsAll = orbitalsAB.QMAtoms();
  bool coordinates_agree = true;
  for (Index i = 0; i < atomsAll.size(); i++) {
    const QMAtom& dimer = atomsAll[i];
    const QMAtom* monomer = nullptr;

    if (i < atomsA.size()) {
      monomer = &atomsA[i];
    } else if (i < atomsB.size() + atomsA.size()) {
      monomer = &atomsB[i - atomsA.size()];
    } else {
      // Linker
      XTP_LOG(Log::error, *_pLog)
          << (boost::format(
                  "Neither Monomer A nor Monomer B contains "
                  "atom %s on line %u. Hence, this atom is part of a linker.") %
              dimer.getElement() % (i + 1))
                 .str()
          << std::flush;
      continue;
    }

    if (!monomer->getPos().isApprox(dimer.getPos(), 0.001)) {
      coordinates_agree = false;
    }

    if (monomer->getElement() != dimer.getElement()) {
      throw std::runtime_error(
          "Atom types do not agree in dimer and monomers\n");
    }
  }

  if (!coordinates_agree) {
    XTP_LOG(Log::error, *_pLog)
        << "======WARNING=======\n Coordinates of monomer "
           "and dimer atoms do not agree"
        << std::flush;
  }
}

}  // namespace xtp
}  // namespace votca
