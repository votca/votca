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
#ifndef VOTCA_XTP_DFTCOUPLING_H
#define VOTCA_XTP_DFTCOUPLING_H

#include <votca/xtp/couplingbase.h>

namespace votca {
namespace xtp {

/**
 * \brief Evaluates electronic coupling elements
 *
 * B. Baumeier, J. Kirkpatrick, D. Andrienko,
 * Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
 *
 */

class DFTcoupling : public CouplingBase {
 public:
  std::string Identify() const { return "dftcoupling"; }

  void CalculateCouplings(const Orbitals& orbitalsA, const Orbitals& orbitalsB,
                          const Orbitals& orbitalsAB) override;

  void Initialize(tools::Property&) override;

  void Addoutput(tools::Property& type_summary, const Orbitals& orbitalsA,
                 const Orbitals& orbitalsB) const override;

 private:
  void WriteToProperty(tools::Property& type_summary, const Orbitals& orbitalsA,
                       const Orbitals& orbitalsB, Index a, Index b) const;
  double getCouplingElement(Index levelA, Index levelB,
                            const Orbitals& orbitalsA,
                            const Orbitals& orbitalsB) const;

  std::pair<int, Index> DetermineRangeOfStates(const Orbitals& orbital,
                                               Index numberofstates) const;

  Eigen::MatrixXd JAB;

  double _degeneracy = 0.0;
  Index _numberofstatesA = 1;
  Index _numberofstatesB = 1;

  std::pair<int, Index> Range_orbA;
  std::pair<int, Index> Range_orbB;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTCOUPLING_H
