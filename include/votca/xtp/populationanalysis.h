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
#ifndef VOTCA_XTP_POPULATIONANALYSIS_H
#define VOTCA_XTP_POPULATIONANALYSIS_H

#include <votca/xtp/aomatrix.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmfragment.h>

/**
 * \brief Takes a list of atoms, and the corresponding density matrix and puts
 * out a table of wavefunction partial charges
 *
 *
 *
 */

namespace votca {
namespace xtp {

template <bool T>
class Populationanalysis {
 public:
  StaticSegment CalcChargeperAtom(const Orbitals& orbitals,
                                  const QMState& state) const;

  void CalcChargeperFragment(std::vector<QMFragment<BSE_Population> >& frags,
                             const Orbitals& orbitals, QMStateType type) const;

 private:
  Eigen::VectorXd CalcNucChargeperAtom(const QMMolecule& mol) const;

  Eigen::VectorXd CalcElecChargeperAtom(const Eigen::MatrixXd& dmat,
                                        AOOverlap& overlap,
                                        const AOBasis& basis) const;
};

using Mulliken = Populationanalysis<false>;
using Lowdin = Populationanalysis<true>;

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_POPULATIONANALYSIS_H
