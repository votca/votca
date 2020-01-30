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
#ifndef VOTCA_XTP_COUPLINGBASE_H
#define VOTCA_XTP_COUPLINGBASE_H

#include <boost/format.hpp>

#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

/**
 * \brief Base Class to derive DFT and BSE coupling from
 *
 */

class CouplingBase {
 public:
  virtual void CalculateCouplings(const Orbitals& orbitalsA,
                                  const Orbitals& orbitalsB,
                                  const Orbitals& orbitalsAB) = 0;

  virtual void Initialize(tools::Property&) = 0;

  virtual void Addoutput(tools::Property& type_summary,
                         const Orbitals& orbitalsA,
                         const Orbitals& orbitalsB) const = 0;

  void setLogger(Logger* pLog) { _pLog = pLog; }

 protected:
  Logger* _pLog;
  void CheckAtomCoordinates(const Orbitals& orbitalsA,
                            const Orbitals& orbitalsB,
                            const Orbitals& orbitalsAB) const;

  Eigen::MatrixXd CalculateOverlapMatrix(const Orbitals& orbitalsAB) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_COUPLINGBASE_H
